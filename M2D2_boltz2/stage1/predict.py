"""Stage 1: Boltz-2 binding-affinity predictions for every drug x protein pair.

Drugs default to every agent named in the interaction table (all three agent
columns, so growth media and third agents are included); proteins to the whole
E. coli panel. One protein MSA is fetched once and reused for every drug.

    # everything, one GPU
    python predict.py --fetch-msa

    # split the panel across N array tasks (any scheduler), one GPU each
    python predict.py --shard $I --nshards 24

    # add agents later, against the same panel and cached MSAs
    python predict.py --agents GLU GLY VER

A background thread fetches MSAs while the GPU works on proteins already
fetched, so a cold-cache run does not wait for all of them up front. Proteins
are handed to Boltz in batches of --batch-size as their MSAs land.

Re-running skips any pair whose prediction already exists.  Prediction settings
(sampling steps, diffusion samples) are the ones used for the paper and are
fixed in BOLTZ_ARGS below.
"""
import argparse
import glob
import os
import queue
import shutil
import socket
import subprocess
import threading
import time
from pathlib import Path

import pandas as pd
import yaml

# Interaction-table abbreviations that differ from drugs.xlsx.
ALIAS = {"VPM": "VER"}

# Boltz-2 defaults are 200 sampling steps for both passes and 5 affinity samples.
BOLTZ_ARGS = ["--diffusion_samples_affinity", "5",
              "--sampling_steps", "100",
              "--sampling_steps_affinity", "100"]

MSA_SERVER = "https://api.colabfold.com"
_DONE = object()


def safe(s):
    """Filename-safe form of an ID, as Boltz will echo it back in output paths."""
    return "".join(x for x in str(s).strip() if x.isalnum() or x in "._-")


def select_drugs(drugs, interactions, agents):
    """Rows of drugs.xlsx to predict: an explicit --agents list, or every agent
    appearing anywhere in the interaction table."""
    if agents:
        want = set(agents)
    else:
        t = pd.read_excel(interactions, sheet_name="allInteractions")
        want = set(pd.concat([t[c] for c in ("drug_1", "drug_2", "drug_3")]).dropna())
    want = {ALIAS.get(a, a) for a in want}
    missing = sorted(want - set(drugs["abbrev"]))
    if missing:
        raise SystemExit(f"no row in drugs.xlsx for: {', '.join(missing)}")
    sel = drugs[drugs["abbrev"].isin(want)]
    if sel["SMILES"].isna().any():
        raise SystemExit(f"no SMILES for: {', '.join(sel.loc[sel.SMILES.isna(), 'abbrev'])}")
    return sel


def fetch_msa(pid, seq, msa_dir, tries=5):
    """Fetch one protein's MSA from the ColabFold server, with retry."""
    from boltz.main import compute_msa
    for attempt in range(tries):
        try:
            compute_msa(data={f"{pid}_A": seq}, target_id=pid, msa_dir=msa_dir,
                        msa_server_url=MSA_SERVER, msa_pairing_strategy="greedy")
            time.sleep(1)  # be polite to the API
            return True
        except Exception as e:  # noqa: BLE001 - network errors of every kind
            print(f"  [msa] {pid} attempt {attempt + 1}/{tries} failed: {e}", flush=True)
            time.sleep(20)
    return False


def msa_worker(proteins, msa_dir, fetch, out_queue):
    """Ready each protein's MSA in panel order, ahead of the GPU consuming them."""
    for pid, seq in proteins:
        ready = (msa_dir / f"{pid}_A.csv").exists() or (fetch and fetch_msa(pid, seq, msa_dir))
        out_queue.put((pid, seq, ready))
    out_queue.put(_DONE)


def done_pairs(pred_dir):
    """(drug, protein) pairs that already have a prediction JSON under pred_dir."""
    done = set()
    for f in glob.glob(os.path.join(pred_dir, "**", "affinity_affinity_*.json"), recursive=True):
        stem = os.path.basename(f)[len("affinity_affinity_"):-len(".json")]
        drug, _, prot = stem.rpartition("_")
        done.add((drug, prot))
    return done


def run_boltz(yaml_dir, out_dir, devices):
    cmd = ["boltz", "predict", str(yaml_dir), "--out_dir", str(out_dir),
           "--devices", str(devices), *BOLTZ_ARGS]
    env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
    env["PYTHONWARNINGS"] = "ignore"
    env["TORCH_FLOAT32_MATMUL_PRECISION"] = "medium"  # TF32 on Ampere GPUs
    if devices > 1:
        # Pin DDP to a port we know is free: concurrent shards on one node would
        # otherwise race for the same default and fail with "address already in use".
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.bind(("", 0))
            env["MASTER_PORT"] = str(s.getsockname()[1])
        env.update(MASTER_ADDR="127.0.0.1", NCCL_SOCKET_IFNAME="lo,eth0,eth1,en0",
                   GLOO_SOCKET_IFNAME="lo,eth0,eth1,en0", NCCL_IGNORE_DISABLED_MAC_VLAN="1")
    print(f"  $ {' '.join(cmd)}", flush=True)
    return subprocess.run(cmd, env=env).returncode


def write_yamls(pid, seq, drugs, msa, batch_dir, pred_dir, done):
    """One YAML per drug against this protein; returns how many were written."""
    written = 0
    for _, d in drugs.iterrows():
        did = safe(d["Full name"])
        if (did, pid) in done:
            continue
        spec = {"version": 1,
                "sequences": [{"protein": {"id": "A", "sequence": seq, "msa": str(msa)}},
                              {"ligand": {"id": "B", "smiles": str(d["SMILES"]).strip()}}],
                "properties": [{"affinity": {"binder": "B"}}]}
        with open(batch_dir / f"affinity_{did}_{pid}.yaml", "w") as f:
            yaml.dump(spec, f, sort_keys=False)
        # Boltz skips any record whose output directory exists, even when a
        # previous run crashed before writing the JSON; clear it so the pair is
        # genuinely re-predicted.
        for stale in pred_dir.glob(f"boltz_results_*/predictions/affinity_{did}_{pid}"):
            shutil.rmtree(stale)
        written += 1
    return written


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--drugs", default="drugs.xlsx", help="columns: Full name, abbrev, SMILES")
    ap.add_argument("--proteins", default="ecoli_4087_AAseq.xlsx", help="columns: Protein, Sequence")
    ap.add_argument("--interactions", default="drugInteractions_weights.xlsx",
                    help="defines the drug set when --agents is not given")
    ap.add_argument("--agents", nargs="+", help="abbreviations to predict instead")
    ap.add_argument("--msa-dir", default="msa_cache")
    ap.add_argument("--fetch-msa", action="store_true", help="fetch missing MSAs from ColabFold")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--out", default="boltz_out")
    ap.add_argument("--devices", type=int, default=1)
    ap.add_argument("--batch-size", type=int, default=8, metavar="N",
                    help="proteins per Boltz call; 1 predicts each as its MSA lands, "
                         "0 waits and predicts everything in one call")
    ap.add_argument("--queue-depth", type=int, default=16, metavar="N",
                    help="how far MSA fetching may run ahead of prediction")
    ap.add_argument("--yaml-only", action="store_true", help="write inputs, do not run Boltz")
    a = ap.parse_args()

    drugs = select_drugs(pd.read_excel(a.drugs), a.interactions, a.agents)
    prot = pd.read_excel(a.proteins).drop_duplicates(subset=["Protein", "Sequence"])
    prot = prot.iloc[a.shard::a.nshards]
    print(f"{len(drugs)} drugs x {len(prot)} proteins (shard {a.shard}/{a.nshards})", flush=True)

    msa_dir = Path(a.msa_dir).resolve()
    msa_dir.mkdir(parents=True, exist_ok=True)
    out = Path(a.out)
    pred_dir = out / "pred" / f"shard_{a.shard}"
    # Staged inputs are exactly this run's to-do list -- resume is driven by the
    # prediction JSONs, not by these -- so clear anything an earlier run left,
    # otherwise batches accumulate and already-finished pairs get re-predicted.
    yaml_root = out / "yamls" / f"shard_{a.shard}"
    shutil.rmtree(yaml_root, ignore_errors=True)
    yaml_root.mkdir(parents=True)
    done = done_pairs(str(out / "pred"))

    # Fetch MSAs on a background thread so the GPU is not idle waiting on the
    # server; the bounded queue stops it running arbitrarily far ahead.
    todo = [(safe(p["Protein"]), str(p["Sequence"]).strip()) for _, p in prot.iterrows()]
    ready = queue.Queue(maxsize=max(1, a.queue_depth))
    threading.Thread(target=msa_worker, args=(todo, msa_dir, a.fetch_msa, ready),
                     daemon=True).start()

    written = no_msa = skipped = 0
    batch_proteins, batch_pairs, batch_no = 0, 0, 0
    failed = []

    def batch_dir():
        """Staging directory for the current batch, created on first use."""
        d = yaml_root / f"batch_{batch_no:04d}"
        if not d.exists():
            d.mkdir(parents=True)
        return d

    def flush():
        """Predict everything staged in the current batch, then start a new one."""
        nonlocal batch_no, batch_proteins, batch_pairs
        if batch_pairs and not a.yaml_only:
            print(f"[batch {batch_no}] {batch_pairs} pairs over {batch_proteins} proteins",
                  flush=True)
            if run_boltz(batch_dir(), pred_dir, a.devices) != 0:
                failed.append(batch_no)
                print(f"[batch {batch_no}] boltz exited non-zero; continuing", flush=True)
        batch_no += 1
        batch_proteins = batch_pairs = 0

    while True:
        item = ready.get()
        if item is _DONE:
            break
        pid, seq, ok = item
        if not ok:
            no_msa += 1
            continue
        n = write_yamls(pid, seq, drugs, msa_dir / f"{pid}_A.csv", batch_dir(), pred_dir, done)
        skipped += len(drugs) - n
        if n == 0:
            continue
        written += n
        batch_pairs += n
        batch_proteins += 1
        if a.batch_size and batch_proteins >= a.batch_size:
            flush()
    flush()

    print(f"\n{written} pairs predicted, {skipped} already done, "
          f"{no_msa} proteins without an MSA" + ("" if a.fetch_msa else " (pass --fetch-msa)"))
    if failed:
        raise SystemExit(f"{len(failed)} batch(es) failed: {failed}")


if __name__ == "__main__":
    main()

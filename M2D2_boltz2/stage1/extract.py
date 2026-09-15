"""Stage 1: collect Boltz-2 affinity JSONs into drug x protein matrices.

Finds every prediction under --pred regardless of how the run was sharded, and
writes one CSV per metric Boltz reports (affinity_probability_binary,
affinity_pred_value, and the two per-head variants of each).  Rows are drugs,
keyed by full name; columns follow the order of the protein panel file.

    python extract.py --pred boltz_out/pred
    python extract.py --pred boltz_out/pred --drop-incomplete-proteins   # uniform panel
    python extract.py --pred new_run/pred --merge-into results/affinity_probability_binary.csv

A protein that failed for any drug (typically an out-of-memory on the longest
sequences) leaves a blank cell; --drop-incomplete-proteins removes such columns
so every drug shares the same panel, which is what stage 2 expects.
"""
import argparse
import glob
import json
import os
from pathlib import Path

import pandas as pd


def safe(s):
    return "".join(x for x in str(s).strip() if x.isalnum() or x in "._-")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pred", required=True, help="directory searched recursively for prediction JSONs")
    ap.add_argument("--drugs", default="drugs.xlsx")
    ap.add_argument("--proteins", default="ecoli_4087_AAseq.xlsx")
    ap.add_argument("--results", default="results")
    ap.add_argument("--merge-into", metavar="CSV",
                    help="existing matrix: its rows are kept, rows for drugs found here are replaced")
    ap.add_argument("--drop-incomplete-proteins", action="store_true")
    a = ap.parse_args()

    files = glob.glob(os.path.join(a.pred, "**", "affinity_affinity_*.json"), recursive=True)
    print(f"{len(files)} prediction files under {a.pred}")

    # Boltz echoes the YAML stem, affinity_{drug}_{protein}; protein IDs carry no underscore.
    name = {safe(n): n for n in pd.read_excel(a.drugs)["Full name"]}
    panel = pd.read_excel(a.proteins).drop_duplicates(subset=["Protein", "Sequence"])
    panel = [safe(p) for p in panel["Protein"]]

    vals, bad = {}, 0
    for f in files:
        stem = os.path.basename(f)[len("affinity_affinity_"):-len(".json")]
        drug, _, prot = stem.rpartition("_")
        try:
            res = json.load(open(f))
        except (OSError, ValueError):
            bad += 1
            continue
        for k, v in res.items():
            vals.setdefault(k, {}).setdefault(name.get(drug, drug), {})[prot] = v
    if bad:
        print(f"  {bad} files unreadable")
    if not vals:
        raise SystemExit("no predictions found")

    key = "affinity_probability_binary"
    drugs = sorted(vals[key])
    for d in drugs:
        n = sum(p in vals[key][d] for p in panel)
        print(f"  {d:<22} {n}/{len(panel)} proteins")

    incomplete = [p for p in panel if any(p not in vals[key][d] for d in drugs)]
    if incomplete:
        shown = " ".join(incomplete[:20]) + (" ..." if len(incomplete) > 20 else "")
        print(f"\n{len(incomplete)} proteins missing for at least one drug: {shown}")
    cols = [p for p in panel if p not in incomplete] if a.drop_incomplete_proteins else panel

    os.makedirs(a.results, exist_ok=True)
    for k, by_drug in vals.items():
        m = pd.DataFrame.from_dict(by_drug, orient="index").reindex(index=drugs, columns=cols)
        m.index.name = "Drug"
        if a.merge_into:
            base = pd.read_csv(Path(a.merge_into).with_name(f"{k}.csv"), index_col="Drug")
            m = pd.concat([base.drop(index=m.index, errors="ignore"), m.reindex(columns=base.columns)])
        m.to_csv(os.path.join(a.results, f"{k}.csv"))
    print(f"\nwrote {len(vals)} matrices to {a.results}/: {len(drugs)} drugs x {len(cols)} proteins")


if __name__ == "__main__":
    main()

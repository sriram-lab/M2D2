"""Stage 2: train and evaluate DDI models on the Boltz-2 feature basis.

Two experiment groups, run together by default:

  representations  Five ways of turning the affinity matrix into a combination
                   feature -- unbinarized, a global 0.5 cutoff, and a per-drug
                   80th-percentile cutoff, each with sigma alone or sigma and
                   delta -- crossed with Random Forest and XGBoost, under
                   5-fold cross-validation over combinations and under
                   leave-one-drug-out. Establishes which representation and
                   which regressor to carry forward.

  published        The best representation (sigma + delta at the per-drug
                   cutoff) under the two protocols M2D2 itself reports, so the
                   Boltz-2 basis can be read against M2D2's published numbers:
                   50 iterations of a 70/30 random holdout with a 100-tree
                   forest at MATLAB TreeBagger settings, and leave-one-drug-out.

    python evaluate.py --matrix results/affinity_probability_binary.csv
    python evaluate.py --experiment published

Writes one CSV per experiment to --results, plus the per-drug leave-one-drug-out
breakdown behind the representation means.
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.model_selection import KFold, train_test_split

from utils import (FIXED_PARAMS, TREEBAGGER, build_model, design_matrix,
                   eval_metrics, ldo_splits, load_combinations, load_drug_classes,
                   load_profiles)

# name, profile key, whether to append the delta score
REPRESENTATIONS = [
    ("sigma (no binarization)",             "raw",   False),
    ("sigma (fixed cutoff 0.5)",            "fixed", False),
    ("sigma + delta (fixed cutoff 0.5)",    "fixed", True),
    ("sigma (per-drug percentile)",         "pct",   False),
    ("sigma + delta (per-drug percentile)", "pct",   True),
]
BEST = ("pct", True)


def run_representations(lists, y, profiles, classes, seed, models):
    """Every representation x regressor, under 5-fold CV and leave-one-drug-out."""
    summary, per_drug = [], []
    for name, key, with_delta in REPRESENTATIONS:
        X = design_matrix(lists, profiles[key], with_delta)
        for model in models:
            folds = []
            for tr, te in KFold(5, shuffle=True, random_state=seed).split(X):
                p = build_model(model, **FIXED_PARAMS[model]).fit(X[tr], y[tr]).predict(X[te])
                folds.append(eval_metrics(y[te], p))

            recs = []
            for tr, te, drug, cls in ldo_splits(lists, classes):
                p = build_model(model, **FIXED_PARAMS[model]).fit(X[tr], y[tr]).predict(X[te])
                recs.append({"representation": name, "model": model, "drug": drug,
                             "class": cls, "n": int(te.size), **eval_metrics(y[te], p)})
            ldo = pd.DataFrame(recs)
            per_drug.append(ldo)

            summary.append({
                "representation": name, "model": model, "dim": X.shape[1], "n_combinations": len(y),
                "cv_pearson": np.mean([f["Pearson"] for f in folds]),
                "cv_spearman": np.mean([f["Spearman"] for f in folds]),
                "cv_kendall": np.mean([f["Kendall"] for f in folds]),
                "n_drugs": len(ldo),
                "ldo_pearson": ldo.Pearson.mean(), "ldo_pearson_sd": ldo.Pearson.std(),
                "ldo_spearman": ldo.Spearman.mean(), "ldo_spearman_sd": ldo.Spearman.std(),
                "ldo_kendall": ldo.Kendall.mean(), "ldo_kendall_sd": ldo.Kendall.std()})
            s = summary[-1]
            print(f"  {name:38s} {model:3s} 5-fold r={s['cv_pearson']:.4f} "
                  f"rho={s['cv_spearman']:.4f} | LDO r={s['ldo_pearson']:.4f} "
                  f"rho={s['ldo_spearman']:.4f} ({len(ldo)} drugs)", flush=True)
    return pd.DataFrame(summary), pd.concat(per_drug, ignore_index=True)


def run_published(lists, y, profiles, classes, n_iter, models):
    """M2D2's own two protocols, on the best representation."""
    key, with_delta = BEST
    X = design_matrix(lists, profiles[key], with_delta)
    rows = []

    r, rho = [], []
    for i in range(n_iter):
        tr, te = train_test_split(np.arange(len(y)), test_size=0.30, random_state=i)
        p = build_model("rf", random_state=i, **TREEBAGGER).fit(X[tr], y[tr]).predict(X[te])
        r.append(pearsonr(y[te], p)[0])
        rho.append(spearmanr(y[te], p)[0])
    for metric, v in [("Pearson r", r), ("Spearman rho", rho)]:
        rows.append({"experiment": f"70/30 holdout x{n_iter}, RF(100)", "metric": metric,
                     "value": np.mean(v), "sd": np.std(v, ddof=1),
                     "se": np.std(v, ddof=1) / np.sqrt(len(v)), "n_units": len(v),
                     "spread_over": "random 70/30 splits"})
    print(f"  70/30 holdout x{n_iter}  r={np.mean(r):.4f}+-{np.std(r, ddof=1):.4f}  "
          f"rho={np.mean(rho):.4f}", flush=True)

    for model in models:
        label = {"rf": "Random Forest", "xgb": "XGBoost"}[model]
        r, rho = [], []
        for tr, te, _, _ in ldo_splits(lists, classes):
            p = build_model(model, **FIXED_PARAMS[model]).fit(X[tr], y[tr]).predict(X[te])
            r.append(pearsonr(y[te], p)[0])
            rho.append(spearmanr(y[te], p)[0])
        for metric, v in [("Pearson r", r), ("Spearman rho", rho)]:
            rows.append({"experiment": f"Leave-one-drug-out, {label}", "metric": metric,
                         "value": np.mean(v), "sd": np.std(v, ddof=1),
                         "se": np.std(v, ddof=1) / np.sqrt(len(v)), "n_units": len(v),
                         "spread_over": "held-out drugs"})
        print(f"  LDO {label:14s} r={np.mean(r):.4f}+-{np.std(r, ddof=1):.4f}  "
              f"rho={np.mean(rho):.4f} ({len(r)} drugs)", flush=True)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--matrix", default="results/affinity_probability_binary.csv",
                    help="stage-1 output: drugs x proteins")
    ap.add_argument("--drugs", default="drugs.xlsx")
    ap.add_argument("--interactions", default="drugInteractions_weights.xlsx")
    ap.add_argument("--classes", default="leaveOneDrugOut_w75_n25_v2.xlsx")
    ap.add_argument("--experiment", choices=["representations", "published", "all"], default="all")
    ap.add_argument("--percentile", type=float, default=80.0, help="per-drug binarization cutoff")
    ap.add_argument("--iterations", type=int, default=50, help="70/30 holdout repeats")
    ap.add_argument("--seed", type=int, default=0, help="seed for the 5-fold split")
    ap.add_argument("--models", nargs="+", choices=["rf", "xgb"], default=["rf", "xgb"],
                    help="regressors to evaluate; xgb needs xgboost installed")
    ap.add_argument("--results", default="results")
    a = ap.parse_args()

    profiles = load_profiles(a.matrix, a.drugs, a.percentile)
    lists, y = load_combinations(a.interactions, profiles["pct"])
    classes = load_drug_classes(a.classes)
    os.makedirs(a.results, exist_ok=True)

    if a.experiment in ("representations", "all"):
        print("\nrepresentations:")
        summary, per_drug = run_representations(lists, y, profiles, classes, a.seed, a.models)
        summary.to_csv(f"{a.results}/representations.csv", index=False)
        per_drug.to_csv(f"{a.results}/representations_ldo_per_drug.csv", index=False)

    if a.experiment in ("published", "all"):
        print("\npublished protocols:")
        run_published(lists, y, profiles, classes, a.iterations, a.models).to_csv(
            f"{a.results}/published_protocols.csv", index=False)

    print(f"\nwrote results to {a.results}/")


if __name__ == "__main__":
    main()

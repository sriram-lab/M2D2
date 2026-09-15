"""Stage 2: shared loading, feature construction, models and metrics.

The stage-1 affinity matrix is turned into M2D2's sigma and delta scores, the
interaction table into one label per drug combination, and the two together into
the design matrix every experiment in evaluate.py trains on.
"""
import numpy as np
import pandas as pd
from scipy.stats import kendalltau, pearsonr, spearmanr
from sklearn.ensemble import RandomForestRegressor

AGENTS = ["drug_1", "drug_2", "drug_3"]
UNCLASSIFIED = "(unclassified)"

# drugs.xlsx abbreviates verapamil VER; the interaction table calls it VPM.
ALIAS = {"VER": "VPM"}

# Hyperparameters held fixed across every representation and protocol.
FIXED_PARAMS = {
    "rf": {"n_estimators": 400, "max_features": "sqrt", "min_samples_leaf": 1, "max_depth": None},
    "xgb": {"n_estimators": 400, "max_depth": 6, "learning_rate": 0.05, "subsample": 0.8},
}

# MATLAB's TreeBagger regression defaults, which M2D2 trains with; sklearn's
# RandomForestRegressor defaults differ on both counts.
TREEBAGGER = {"n_estimators": 100, "max_features": 1.0 / 3.0, "min_samples_leaf": 5}


def build_model(model_name, device="cpu", random_state=0, **params):
    """A regressor with everything but `params` pinned. Seeded, but note that
    n_jobs=-1 accumulates in nondeterministic order, so repeated runs can differ
    by ~1e-16 per prediction and ~1e-4 in a mean correlation. Pass n_jobs=1 for
    bit-exact repeats.

    The holdout protocol reseeds per iteration, so random_state is overridable.
    """
    if model_name == "rf":
        return RandomForestRegressor(random_state=random_state, n_jobs=-1, **params)
    if model_name == "xgb":
        from xgboost import XGBRegressor  # optional: only the xgb path needs it

        return XGBRegressor(random_state=random_state, tree_method="hist", device=device,
                            n_jobs=8 if device == "cpu" else None, **params)
    raise ValueError(model_name)


def _safe(fn, y_true, y_pred):
    """Correlation that returns 0 rather than NaN when a fold has no variance."""
    stat, p = fn(y_true, y_pred)
    return (0.0, p) if np.isnan(stat) else (float(stat), p)


def eval_metrics(y_true, y_pred):
    if len(y_true) > 1:
        pr, pr_p = _safe(pearsonr, y_true, y_pred)
        sp, sp_p = _safe(spearmanr, y_true, y_pred)
        kt, kt_p = _safe(kendalltau, y_true, y_pred)
    else:
        pr, pr_p = sp, sp_p = kt, kt_p = 0.0, 1.0
    return {"Pearson": pr, "Pearson_pval": pr_p, "Spearman": sp, "Spearman_pval": sp_p,
            "Kendall": kt, "Kendall_pval": kt_p, "Validation size": len(y_true)}


# --------------------------------------------------------------------- features

def load_profiles(matrix, drugs_path, pct=80.0):
    """Read the stage-1 matrix and return the three binarization variants.

    Rows are keyed by the interaction table's abbreviation. Protein columns with
    any missing value are dropped, so every drug shares one panel.
    """
    n2a = dict(zip(*pd.read_excel(drugs_path)[["Full name", "abbrev"]].values.T))
    aff = pd.read_csv(matrix)
    aff["ab"] = aff["Drug"].map(n2a).map(lambda a: ALIAS.get(a, a))
    aff = aff.dropna(subset=["ab"]).dropna(axis=1)
    cols = [c for c in aff.columns if c not in ("Drug", "ab")]
    m = aff.set_index("ab")[cols].astype(float)
    print(f"panel: {m.shape[0]} drugs x {m.shape[1]} proteins")
    return {
        "raw": {d: v.to_numpy() for d, v in m.iterrows()},
        "fixed": {d: v.to_numpy() for d, v in (m > 0.5).astype(float).iterrows()},
        "pct": {d: v.to_numpy() for d, v in
                m.gt(m.quantile(pct / 100.0, axis=1), axis=0).astype(float).iterrows()},
    }


def joint(agents, prof, with_delta):
    """M2D2's sigma (union) and delta (unique) scores for a combination.

    Both are symmetric in the agents and defined for any number of them, so
    pairs and triples produce features of the same length.
    """
    hits = np.sum([prof[a] for a in agents], axis=0)
    sigma = hits * (2.0 / len(agents))
    return np.concatenate([sigma, (hits == 1).astype(float)]) if with_delta else sigma


def design_matrix(lists, prof, with_delta):
    return np.vstack([joint(l, prof, with_delta) for l in lists])


# ----------------------------------------------------------------------- labels

def _dedup(df, cols):
    """One row per unordered combination, keeping the larger score.

    The source table scores 87 pairs under both Bliss and Loewe. Identical agents
    give identical features, so leaving both copies in lets a random split put one
    in training and its twin in test.
    """
    key = df.apply(lambda r: tuple(sorted(r[c] for c in cols)), axis=1)
    g = df.assign(_k=key).groupby("_k", as_index=False)["score"].max()
    for i, c in enumerate(cols):
        g[c] = g["_k"].map(lambda k, i=i: k[i])
    return g[cols + ["score"]]


def load_combinations(interactions, profiles, sheet="avgBliss_avgLoewe"):
    """Pairs and three-agent combinations pooled into one labelled set.

    Returns (lists, y): the agents of each combination, and its synergy score.
    The default sheet holds one row per measurement with Bliss and Loewe scores
    already averaged within each; its pairwise rows are what earlier versions of
    this pipeline kept as a separate processed_drug_interactions.csv, and its
    triples are identical to the allInteractions sheet's.
    """
    t = pd.read_excel(interactions, sheet_name=sheet)
    is_triple = t.drug_3.notna() & (t.drug_3.astype(str).str.strip().str.lower() != "nan")

    pairs = t[~is_triple]
    pairs = pairs[[all(a in profiles for a in (r.drug_1, r.drug_2)) for r in pairs.itertuples()]]
    pairs = _dedup(pairs, AGENTS[:2]).assign(drug_3=np.nan)

    trips = t[is_triple]
    trips = trips[[all(a in profiles for a in (r.drug_1, r.drug_2, r.drug_3))
                   for r in trips.itertuples()]]
    trips = _dedup(trips, AGENTS)

    pooled = pd.concat([pairs, trips], ignore_index=True)
    lists = [[r[c] for c in AGENTS if pd.notna(r[c]) and str(r[c]).strip() != ""]
             for _, r in pooled.iterrows()]
    print(f"combinations: {len(pairs)} pairs + {len(trips)} triples = {len(pooled)}")
    return lists, pooled["score"].to_numpy(float)


def load_drug_classes(path):
    """Held-out drug -> mechanistic class, for the leave-one-drug-out breakdown."""
    t = pd.read_excel(path, sheet_name="ml (2)")[["allDrugs_nplus2", "Unnamed: 8"]].dropna()
    return dict(zip(t["allDrugs_nplus2"], t["Unnamed: 8"]))


def ldo_splits(lists, drug_classes, min_test_size=25):
    """Leave-one-drug-out splits, over drugs carrying a mechanistic class.

    M2D2 requires a held-out drug to appear in at least 25 interactions so the
    test fold supports a stable correlation.
    """
    n = len(lists)
    for d in sorted({a for l in lists for a in l}):
        if d not in drug_classes:
            continue
        test = np.array([i for i, l in enumerate(lists) if d in l])
        if min_test_size <= test.size < n:
            yield np.setdiff1d(np.arange(n), test), test, d, drug_classes.get(d, UNCLASSIFIED)

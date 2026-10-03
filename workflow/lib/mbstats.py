"""Shared helpers for the IBD cross-platform microbiome workflow.

Everything here is plain pandas/numpy so that it can be unit-tested without the
heavy bioinformatics stack (see tests/test_units.py).
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

REQUIRED_COLUMNS = [
    "sample_id", "specimen_id", "subject_id", "data_type", "batch", "sample_type",
    "week", "diagnosis", "fastq_1", "fastq_2",
]
NUMERIC_COLUMNS = ["week", "age", "bmi", "hbi", "sccai"]
YES = {"yes", "y", "true", "1"}
NO = {"no", "n", "false", "0"}
NA_STRINGS = ["", "NA", "NaN", "nan", "N/A", "n/a", "None", "NULL"]


# --------------------------------------------------------------------------
# I/O
# --------------------------------------------------------------------------
def read_layer(path) -> pd.DataFrame:
    """Read a features x samples table whose first column holds feature names."""
    df = pd.read_csv(path, sep="\t", index_col=0, na_values=NA_STRINGS[1:])
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df.astype(float)


def write_layer(df: pd.DataFrame, path, index_label="feature") -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index_label=index_label, float_format="%.10g")


# --------------------------------------------------------------------------
# Sample sheet
# --------------------------------------------------------------------------
def _norm_yes_no(x):
    if pd.isna(x):
        return np.nan
    s = str(x).strip().lower()
    if s in YES:
        return "Yes"
    if s in NO:
        return "No"
    return str(x)


def load_samples(path, diagnosis_levels=("nonIBD", "CD", "UC")) -> pd.DataFrame:
    """Read + validate the sample sheet. Raises ValueError with an actionable message."""
    df = pd.read_csv(path, sep="\t", dtype=str, na_values=NA_STRINGS, keep_default_na=False)
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing required column(s): {missing}")
    for c in NUMERIC_COLUMNS:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    for c in ("antibiotics", "immunosuppressants"):
        if c in df.columns:
            df[c] = df[c].map(_norm_yes_no)
    if "assemble" in df.columns:
        df["assemble"] = df["assemble"].map(lambda v: str(v).strip().lower() in YES)
    else:
        df["assemble"] = False

    if df["sample_id"].isna().any() or df["sample_id"].duplicated().any():
        dup = df.loc[df["sample_id"].duplicated(), "sample_id"].tolist()
        raise ValueError(f"sample_id must be unique and non-empty; duplicates: {dup[:5]}")
    if not df["sample_id"].str.fullmatch(r"[A-Za-z0-9_.\-]+").all():
        raise ValueError("sample_id may only contain letters, digits, '_', '.' and '-'")
    bad_type = set(df["data_type"]) - {"16S", "MGX"}
    if bad_type:
        raise ValueError(f"data_type must be '16S' or 'MGX'; found {bad_type}")
    bad_st = set(df["sample_type"]) - {"sample", "negative_control"}
    if bad_st:
        raise ValueError(f"sample_type must be 'sample' or 'negative_control'; found {bad_st}")

    real = df[df["sample_type"] == "sample"]
    for c in ("specimen_id", "subject_id", "diagnosis", "week", "batch"):
        if real[c].isna().any():
            raise ValueError(f"column '{c}' has missing values for sample_type == 'sample'")
    bad_dx = set(real["diagnosis"]) - set(diagnosis_levels)
    if bad_dx:
        raise ValueError(f"unknown diagnosis labels {bad_dx}; allowed: {list(diagnosis_levels)}")
    dup = real.duplicated(["specimen_id", "data_type"])
    if dup.any():
        raise ValueError("a specimen may have at most one sample per data_type; offending: "
                         f"{real.loc[dup, 'specimen_id'].tolist()[:5]}")
    # specimen-level fields must agree across platforms
    for c in ("subject_id", "diagnosis", "week"):
        n = real.groupby("specimen_id")[c].nunique(dropna=False)
        if (n > 1).any():
            raise ValueError(f"'{c}' differs between samples of the same specimen: {n[n > 1].index.tolist()[:5]}")
    # one diagnosis per subject
    n = real.groupby("subject_id")["diagnosis"].nunique()
    if (n > 1).any():
        raise ValueError(f"subjects with >1 diagnosis: {n[n > 1].index.tolist()[:5]}")
    if df[["fastq_1", "fastq_2"]].isna().any().any():
        raise ValueError("fastq_1/fastq_2 required for every row (paired-end data only)")
    return df.reset_index(drop=True)


# --------------------------------------------------------------------------
# Clinical derivations: activity, pre-flare label, flare events
# --------------------------------------------------------------------------
def derive_clinical(samples: pd.DataFrame, design: dict):
    """Return (sample_meta, specimens, flare_events).

    * specimens: one row per biological specimen with activity/flare labels.
    * flare_events: visits aligned to the first active visit after an inactive one.
    * sample_meta: per-sample metadata (+ specimen-level labels, is_baseline, visit_index).
    """
    ref = design.get("reference_level", "nonIBD")
    act = design["activity"]
    gap = float(design["flare"]["max_gap_weeks"])
    window = list(range(design["flare"]["window"][0], design["flare"]["window"][1] + 1))

    real = samples[samples["sample_type"] == "sample"].copy()
    spec_cols = [c for c in real.columns
                 if c not in ("sample_id", "data_type", "batch", "fastq_1", "fastq_2", "assemble", "sample_type")]
    specimens = (real.sort_values("data_type")
                 .groupby("specimen_id", as_index=False)[spec_cols]
                 .first())
    specimens = specimens.drop(columns=[c for c in specimens.columns if c == "specimen_id.1"], errors="ignore")
    plat = real.pivot(index="specimen_id", columns="data_type", values="sample_id")
    plat.columns = [f"sample_{c.lower()}" for c in plat.columns]
    specimens = specimens.merge(plat, left_on="specimen_id", right_index=True, how="left")

    def activity(row):
        cfg = act.get(row["diagnosis"])
        if cfg is None or cfg["column"] not in row.index:
            return np.nan
        v = row[cfg["column"]]
        return np.nan if pd.isna(v) else float(v >= cfg["threshold"])

    specimens["active"] = specimens.apply(activity, axis=1)
    specimens = specimens.sort_values(["subject_id", "week", "specimen_id"]).reset_index(drop=True)
    specimens["visit_index"] = specimens.groupby("subject_id").cumcount() + 1

    g = specimens.groupby("subject_id")
    specimens["next_week"] = g["week"].shift(-1)
    specimens["next_active"] = g["active"].shift(-1)
    specimens["prev_week"] = g["week"].shift(1)
    specimens["prev_active"] = g["active"].shift(1)
    linked = (specimens["next_week"] - specimens["week"]) <= gap
    flare = pd.Series(np.nan, index=specimens.index)
    elig = (specimens["diagnosis"] != ref) & (specimens["active"] == 0) & linked & specimens["next_active"].notna()
    flare[elig] = (specimens.loc[elig, "next_active"] == 1).astype(float)
    specimens["flare_next"] = flare

    # flare events
    rows, eid = [], 0
    for subj, d in specimens.groupby("subject_id"):
        d = d.reset_index(drop=True)
        if (d["diagnosis"] == ref).all():
            continue
        for i in range(1, len(d)):
            onset = (d.loc[i - 1, "active"] == 0 and d.loc[i, "active"] == 1
                     and (d.loc[i, "week"] - d.loc[i - 1, "week"]) <= gap)
            if not onset:
                continue
            eid += 1
            for rel in window:
                j = i + rel
                if 0 <= j < len(d):
                    rows.append({"event_id": f"{subj}_E{eid}", "subject_id": subj,
                                 "specimen_id": d.loc[j, "specimen_id"], "rel_visit": rel,
                                 "week": d.loc[j, "week"], "active": d.loc[j, "active"]})
    events = pd.DataFrame(rows, columns=["event_id", "subject_id", "specimen_id", "rel_visit", "week", "active"])

    keep = ["specimen_id", "visit_index", "active", "flare_next"]
    sample_meta = real.merge(specimens[keep], on="specimen_id", how="left")
    first = sample_meta.groupby(["subject_id", "data_type"])["week"].transform("min")
    sample_meta["is_baseline"] = sample_meta["week"] == first
    # break ties (identical weeks) deterministically
    sample_meta = sample_meta.sort_values(["subject_id", "data_type", "week", "sample_id"])
    dup = sample_meta.duplicated(["subject_id", "data_type", "is_baseline"]) & sample_meta["is_baseline"]
    sample_meta.loc[dup, "is_baseline"] = False
    return sample_meta.reset_index(drop=True), specimens, events


# --------------------------------------------------------------------------
# Compositional helpers
# --------------------------------------------------------------------------
def to_proportions(df: pd.DataFrame) -> pd.DataFrame:
    s = df.sum(axis=0)
    return df.div(s.where(s > 0), axis=1).fillna(0.0)


def clr(df: pd.DataFrame, kind: str = "counts") -> pd.DataFrame:
    """CLR transform of a features x samples table.

    counts   : add a pseudocount of 0.5 to every cell.
    abundance: close to proportions, add half the smallest non-zero proportion, re-close.
    """
    x = df.to_numpy(dtype=float)
    if kind == "counts":
        x = x + 0.5
    else:
        p = x / np.where(x.sum(0) > 0, x.sum(0), 1.0)
        pos = p[p > 0]
        pc = pos.min() / 2 if pos.size else 1e-6
        x = p + pc
    x = x / x.sum(0, keepdims=True)
    lx = np.log(x)
    return pd.DataFrame(lx - lx.mean(0, keepdims=True), index=df.index, columns=df.columns)


def filter_features(df: pd.DataFrame, min_prev: float, min_mean_prop: float = 0.0) -> pd.DataFrame:
    prev = (df > 0).mean(axis=1)
    keep = prev >= min_prev
    if min_mean_prop > 0:
        keep &= to_proportions(df).mean(axis=1) >= min_mean_prop
    return df.loc[keep]


# --------------------------------------------------------------------------
# Taxon-name harmonisation (SILVA / GTDB / MetaPhlAn)
# --------------------------------------------------------------------------
_UNCLASSIFIED = re.compile(r"^(unclassified|uncultured|unknown|incertae|ggb\d+|sgb\d+|\s*$)", re.I)


def normalize_taxon_name(name) -> str | None:
    """Reduce a genus label to a comparable key; None if unusable.

    'g__Blautia_A' -> 'blautia'; '[Ruminococcus] gnavus group' -> 'ruminococcus';
    'Escherichia-Shigella' -> 'escherichia'; 'Unclassified_Lachnospiraceae' -> None.
    """
    if name is None or (isinstance(name, float) and np.isnan(name)):
        return None
    s = str(name).strip()
    s = re.sub(r"^[a-z]__", "", s)
    s = s.replace("[", "").replace("]", "")
    if _UNCLASSIFIED.match(s) or s.lower().startswith("unclassified"):
        return None
    s = re.split(r"[ _\-/]", s)[0] if re.match(r"^[A-Za-z]+[ _\-/]", s) else s
    s = re.sub(r"_[A-Z]{1,2}$", "", s)
    s = re.sub(r"[^A-Za-z]", "", s).lower()
    return s or None


def collapse_by_key(df: pd.DataFrame, keys: dict) -> pd.DataFrame:
    """Sum rows sharing the same non-null key."""
    k = pd.Series({f: keys.get(f) for f in df.index}).dropna()
    return df.loc[k.index].groupby(k.values).sum()


# --------------------------------------------------------------------------
# SparCC (Friedman & Alm 2012)
# --------------------------------------------------------------------------
def _basis_corr(frac: np.ndarray, th: float = 0.1, xiter: int = 10):
    """One SparCC iteration on a (samples x features) fraction matrix."""
    n, p = frac.shape
    logf = np.log(frac)
    # variation matrix V_ij = var(log(x_i / x_j))
    V = np.empty((p, p))
    for i in range(p):
        V[i] = np.var(logf[:, [i]] - logf, axis=0, ddof=1)
    D = np.ones((p, p)) + np.diag(np.full(p, p - 2.0))
    excl = np.zeros((p, p), dtype=bool)
    np.fill_diagonal(excl, True)
    cor = np.zeros((p, p))
    for it in range(xiter + 1):
        mask = ~excl
        Dm = mask.astype(float)
        np.fill_diagonal(Dm, mask.sum(1))
        t = (V * mask).sum(1)
        w = np.linalg.lstsq(Dm, t, rcond=None)[0]
        w = np.clip(w, 1e-4, None)
        cov = (w[:, None] + w[None, :] - V) / 2.0
        cor = np.clip(cov / np.sqrt(np.outer(w, w)), -1, 1)
        np.fill_diagonal(cor, 1.0)
        if it == xiter:
            break
        c = np.where(excl, 0, np.abs(cor))
        i, j = np.unravel_index(np.argmax(c), c.shape)
        if c[i, j] <= th:
            break
        excl[i, j] = excl[j, i] = True
    return cor


def sparcc(counts: np.ndarray, iterations: int = 20, rng=None, th: float = 0.1, xiter: int = 10):
    """Median-of-Dirichlet-resamples SparCC correlation. counts: samples x features."""
    rng = np.random.default_rng(rng)
    cors = []
    for _ in range(iterations):
        g = rng.gamma(counts + 1.0)
        frac = g / g.sum(1, keepdims=True)
        cors.append(_basis_corr(frac, th=th, xiter=xiter))
    return np.median(np.stack(cors), axis=0)


def sparcc_with_fdr(counts: np.ndarray, iterations=20, permutations=100, seed=1):
    """Return (cor, fdr) using a pooled permutation null (feature-wise shuffles)."""
    rng = np.random.default_rng(seed)
    obs = sparcc(counts, iterations, rng)
    p = obs.shape[0]
    iu = np.triu_indices(p, 1)
    null = []
    for _ in range(permutations):
        perm = np.apply_along_axis(rng.permutation, 0, counts)
        null.append(np.abs(_basis_corr(_dirichlet_frac(perm, rng)))[iu])
    null = np.sort(np.concatenate(null))
    r = np.abs(obs[iu])
    order = np.argsort(-r)
    rs = r[order]
    exp = (len(null) - np.searchsorted(null, rs, side="left")) / permutations
    fdr = np.minimum(1.0, exp / np.arange(1, len(rs) + 1))
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]   # q-value style monotonicity
    q = np.empty_like(fdr)
    q[order] = fdr
    Q = np.ones((p, p))
    Q[iu] = q
    Q = np.minimum(Q, Q.T)
    return obs, Q


def _dirichlet_frac(counts, rng):
    g = rng.gamma(counts + 1.0)
    return g / g.sum(1, keepdims=True)

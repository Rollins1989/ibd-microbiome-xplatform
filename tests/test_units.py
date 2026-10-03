"""Unit tests for workflow/lib/mbstats.py  (run: pytest -q tests)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "lib"))
import mbstats as mb  # noqa: E402

DESIGN = {
    "reference_level": "nonIBD",
    "activity": {"CD": {"column": "hbi", "threshold": 5}, "UC": {"column": "sccai", "threshold": 5}},
    "flare": {"max_gap_weeks": 8, "window": [-1, 1]},
}


def _sheet():
    rows = []
    # subject P1 (CD): inactive, inactive, ACTIVE (flare at visit 3), active
    for i, (wk, hbi) in enumerate([(0, 2), (4, 3), (8, 7), (12, 8)]):
        rows.append(dict(sample_id=f"P1_{i}_MGX", specimen_id=f"P1_{i}", subject_id="P1", data_type="MGX",
                         batch="b", sample_type="sample", week=wk, diagnosis="CD", hbi=hbi, sccai=np.nan,
                         fastq_1="a", fastq_2="b"))
    # subject P2 (UC): inactive, then a gap of 20 weeks, then active (NOT linked)
    for i, (wk, sc) in enumerate([(0, 1), (20, 9)]):
        rows.append(dict(sample_id=f"P2_{i}_MGX", specimen_id=f"P2_{i}", subject_id="P2", data_type="MGX",
                         batch="b", sample_type="sample", week=wk, diagnosis="UC", hbi=np.nan, sccai=sc,
                         fastq_1="a", fastq_2="b"))
    # control
    rows.append(dict(sample_id="C_0_MGX", specimen_id="C_0", subject_id="C", data_type="MGX", batch="b",
                     sample_type="sample", week=0, diagnosis="nonIBD", hbi=np.nan, sccai=np.nan,
                     fastq_1="a", fastq_2="b"))
    return pd.DataFrame(rows)


def test_flare_labels_and_events():
    meta, spec, ev = mb.derive_clinical(_sheet(), DESIGN)
    s = spec.set_index("specimen_id")
    assert s.loc["P1_0", "flare_next"] == 0      # inactive -> inactive
    assert s.loc["P1_1", "flare_next"] == 1      # inactive -> active (pre-flare visit)
    assert np.isnan(s.loc["P1_2", "flare_next"])  # already active: not eligible
    assert np.isnan(s.loc["P2_0", "flare_next"])  # next visit too far away
    assert np.isnan(s.loc["C_0", "flare_next"])   # controls have no activity
    assert list(ev["rel_visit"]) == [-1, 0, 1]
    assert list(ev["specimen_id"]) == ["P1_1", "P1_2", "P1_3"]
    assert meta.groupby("subject_id")["is_baseline"].sum().eq(1).all()


def test_clr_rows_sum_to_zero_and_are_scale_invariant():
    rng = np.random.default_rng(0)
    df = pd.DataFrame(rng.poisson(20, (30, 8)), columns=[f"s{i}" for i in range(8)]).astype(float)
    c = mb.clr(df, "counts")
    assert np.allclose(c.sum(axis=0), 0)
    a1 = mb.clr(df.div(df.sum(axis=0), axis=1), "abundance")
    a2 = mb.clr(df.div(df.sum(axis=0), axis=1) * 1e6, "abundance")
    assert np.allclose(a1.values, a2.values)


@pytest.mark.parametrize("raw,expected", [
    ("g__Blautia_A", "blautia"), ("[Ruminococcus] gnavus group", "ruminococcus"),
    ("Escherichia-Shigella", "escherichia"), ("Unclassified_Lachnospiraceae", None),
    ("g__GGB1234", None), ("Faecalibacterium", "faecalibacterium"), (np.nan, None),
])
def test_normalize_taxon_name(raw, expected):
    assert mb.normalize_taxon_name(raw) == expected


def test_sparcc_recovers_planted_correlation():
    rng = np.random.default_rng(1)
    n, p = 120, 12
    base = rng.normal(size=(n, p))
    base[:, 1] = base[:, 0] * 0.9 + rng.normal(scale=0.3, size=n)   # strong +
    base[:, 3] = -base[:, 2] * 0.9 + rng.normal(scale=0.3, size=n)  # strong -
    abund = np.exp(base)
    comp = abund / abund.sum(1, keepdims=True)
    counts = np.vstack([rng.multinomial(20000, r) for r in comp]).astype(float)
    cor, q = mb.sparcc_with_fdr(counts, iterations=5, permutations=20, seed=3)
    assert cor[0, 1] > 0.5 and q[0, 1] <= 0.05
    assert cor[2, 3] < -0.5 and q[2, 3] <= 0.05
    off = np.abs(cor[np.triu_indices(p, 1)])
    top2 = np.sort(off)[-2:]
    assert top2.min() > 0.5           # the two planted pairs are the strongest signals
    assert np.sort(off)[-3] < 0.4     # everything else is background noise
    nulls = [cor[i, j] for i in range(p) for j in range(i + 1, p) if (i, j) not in [(0, 1), (2, 3)]]
    assert np.median(np.abs(nulls)) < 0.25


def test_sample_sheet_validation(tmp_path):
    ex = Path(__file__).resolve().parents[1] / "config" / "samples.example.tsv"
    df = mb.load_samples(ex)
    assert len(df) == 7 and df["assemble"].sum() == 2
    bad = pd.read_csv(ex, sep="\t", dtype=str)
    bad.loc[0, "diagnosis"] = "Crohns"
    bad.to_csv(tmp_path / "bad.tsv", sep="\t", index=False)
    with pytest.raises(ValueError, match="unknown diagnosis"):
        mb.load_samples(tmp_path / "bad.tsv")

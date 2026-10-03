#!/usr/bin/env python3
"""Assemble a short Markdown summary of the key quantitative results from the result tables."""
import argparse
from pathlib import Path

import pandas as pd


def rd(p, **kw):
    p = Path(p)
    return pd.read_csv(p, sep="\t", **kw) if p.exists() and p.stat().st_size else None


def md(df, floatfmt=3):
    if df is None or df.empty:
        return "_not available_\n"
    d = df.copy()
    for c in d.select_dtypes("float"):
        d[c] = d[c].map(lambda v: f"{v:.{floatfmt}g}" if pd.notna(v) else "NA")
    head = "| " + " | ".join(map(str, d.columns)) + " |\n|" + "---|" * len(d.columns) + "\n"
    return head + "\n".join("| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    o = Path(a.outdir)
    s = ["# Analysis summary (auto-generated)\n"]
    s += ["## Cohort\n", md(rd(o / "design/design_summary.tsv")), md(rd(o / "design/specimens_by_diagnosis_platform.tsv"))]
    s += ["## Retained samples per layer\n", md(rd(o / "layers/layer_summary.tsv")), "Dropped samples:\n", md(rd(o / "layers/dropped_samples.tsv"))]
    pw = rd(o / "diversity/permanova_pairwise.tsv")
    s += ["## Beta diversity (PERMANOVA, one baseline sample per subject)\n", md(pw)]
    al = rd(o / "diversity/alpha_models.tsv")
    s += ["## Alpha diversity (mixed models)\n", md(al)]
    s += ["## Cross-platform concordance\n", md(rd(o / "concordance/genus_beta_concordance.tsv")),
          md(rd(o / "concordance/bland_altman_shannon.tsv")), md(rd(o / "concordance/genus_harmonisation_coverage.tsv"))]
    rows = []
    for f in sorted((o / "da").glob("*.consensus.tsv")):
        d = pd.read_csv(f, sep="\t")
        for (layer, c), g in d.groupby(["layer", "contrast"]):
            rows.append({"layer": layer, "contrast": c, "features_tested": len(g), "consensus_hits": int(g.consensus.sum()),
                         "up": int(((g.consensus) & (g.direction == "up")).sum()), "down": int(((g.consensus) & (g.direction == "down")).sum())})
    s += ["## Differential abundance (consensus across methods)\n", md(pd.DataFrame(rows))]
    ml = rd(o / "ml/ml_summary.tsv")
    if ml is not None:
        best = ml.sort_values("auroc_oof_mean_pred", ascending=False).groupby(["task", "feature_set"]).head(1)
        s += ["## Prediction (subject-grouped nested CV; best model per task x feature set)\n",
              md(best[["task", "feature_set", "model", "auroc_mean", "auroc_sd", "auroc_oof_mean_pred", "ci_low", "ci_high"]].sort_values(["task", "feature_set"]))]
    s += ["## MAGs\n", md(rd(o / "mags/analysis/mag_summary_stats.tsv"))]
    s += ["## Networks\n", md(rd(o / "networks/network_metrics.tsv").drop(columns=["top_hubs"], errors="ignore") if rd(o / "networks/network_metrics.tsv") is not None else None)]
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text("\n".join(s))


if __name__ == "__main__":
    main()

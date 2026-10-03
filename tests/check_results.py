#!/usr/bin/env python3
"""Verify that the downstream workflow recovered the signal planted by tests/simulate_data.py."""
import sys
from pathlib import Path

import pandas as pd

R = Path(sys.argv[1] if len(sys.argv) > 1 else "tests/sim/results")
fails = []


def check(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg)
    if not cond:
        fails.append(msg)


# --- differential abundance (consensus tables)
d = pd.read_csv(R / "da/16s_genus.consensus.tsv", sep="\t")
hit = lambda f, c, dr: ((d.feature.str.contains(f, regex=False)) & (d.contrast == c) & d.consensus & (d.direction == dr)).any()  # noqa: E731
check(hit("Escherichia", "CD_vs_nonIBD", "up"), "16S: Escherichia up in CD")
check(hit("Faecalibacterium", "CD_vs_nonIBD", "down"), "16S: Faecalibacterium down in CD")
check(hit("Roseburia", "CD_vs_nonIBD", "down"), "16S: Roseburia down in CD")
check(hit("Escherichia", "UC_vs_nonIBD", "up"), "16S: Escherichia up in UC")
m = pd.read_csv(R / "da/mgx_species.consensus.tsv", sep="\t")
check(((m.feature == "Escherichia coli") & (m.contrast == "CD_vs_nonIBD") & m.consensus & (m.direction == "up")).any(),
      "MGX: E. coli up in CD")
null_hits = d[(d.consensus) & d.feature.isin(["Bifidobacterium", "Collinsella", "Prevotella", "Alistipes"])]
check(len(null_hits) <= 1, f"16S: <=1 false positive among 4 null genera (got {len(null_hits)})")

# --- concordance
c = pd.read_csv(R / "concordance/genus_beta_concordance.tsv", sep="\t")
check((c.iloc[:, 2] > 0.8).all() and (c.p < 0.01).all(), "16S vs MGX: Procrustes/Mantel strongly concordant")
check(pd.read_csv(R / "concordance/bland_altman_shannon.tsv", sep="\t").spearman.iloc[0] > 0.7, "Shannon 16S vs MGX correlated")

# --- beta diversity
p = pd.read_csv(R / "diversity/permanova_pairwise.tsv", sep="\t")
row = p[(p.platform == "16S") & (p.distance == "aitchison") & (p.contrast == "CD vs nonIBD")].iloc[0]
check(row.p <= 0.05, f"PERMANOVA 16S Aitchison CD vs nonIBD significant (p={row.p})")

# --- decontam
rem = pd.read_csv(R / "16s/decontam_removed_asvs.tsv", sep="\t")
check({"ASV9001", "ASV9002", "ASV9003"} <= set(rem.ASV), "decontam removed the 3 planted contaminant ASVs")
check(not any(f.startswith("Ralstonia") for f in pd.read_csv(R / "layers/16s_genus.tsv", sep="\t").feature),
      "no Ralstonia (contaminant genus) in the 16S genus layer")

# --- machine learning
ml = pd.read_csv(R / "ml/ml_summary.tsv", sep="\t")
best = ml[ml.task == "CD_vs_nonIBD"].auroc_oof_mean_pred.max()
check(best > 0.8, f"ML CD vs nonIBD AUROC > 0.8 (best {best:.2f})")
pf = ml[ml.task == "preflare"].auroc_oof_mean_pred
check(pf.between(0.5, 0.85).all(), f"pre-flare AUROC modest as simulated ({pf.min():.2f}-{pf.max():.2f})")
shap = pd.read_csv(R / "ml/ml_shap_importance.tsv", sep="\t")
top = shap[shap.task == "CD_vs_nonIBD"].head(8).feature.str.cat(sep=" ")
check("Escherichia" in top or "Faecalibacterium" in top, "SHAP top features include planted taxa")

# --- longitudinal
a = pd.read_csv(R / "longitudinal/activity_associations.tsv", sep="\t")
check(((a.feature.str.contains("Escherichia")) & (a.estimate > 0) & (a.q < 0.05)).any(), "activity: Escherichia higher when active")

# --- MAGs
ms = pd.read_csv(R / "mags/analysis/mag_summary_stats.tsv", sep="\t").iloc[0]
check(ms.n_novel_species >= 1 and ms.n_mags == 36, "MAG catalogue annotated, novel species flagged")

# --- figures
figs = ["Fig1_study_design", "Fig2_diversity", "Fig3_platform_concordance", "Fig4_differential_abundance",
        "Fig5_ml_performance", "Fig6_MAG_tree", "Fig7_longitudinal", "Fig8_cooccurrence_networks", "Fig9_strain_persistence"]
for f in figs:
    check((R / "figures" / f"{f}.pdf").stat().st_size > 5000, f"{f}.pdf exists and is non-trivial")

print(f"\n{len(fails)} failed check(s)")
sys.exit(1 if fails else 0)

#!/usr/bin/env python3
"""StrainPhlAn follow-up: within-subject strain persistence vs between-subject diversity.

For every clade tree: patristic tip-to-tip distances; pairs of samples from the same subject
(within) are compared with pairs from different subjects (between). Persistence = within << between.
Within-subject distances are also summarised by diagnosis and time gap.
Samples of the same subject are NOT independent of each other, so p-values compare pooled distance
distributions descriptively (Mann-Whitney) - use them as a sanity check, not as inference.
"""
import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402
from skbio import TreeNode  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
import mbstats as mb  # noqa: E402

COLORS = {"nonIBD": "#0072B2", "CD": "#D55E00", "UC": "#009E73"}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--trees", nargs="*", default=[])
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--figdir", required=True)
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    Path(a.figdir).mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(a.meta, sep="\t", na_values=mb.NA_STRINGS, keep_default_na=False).set_index("sample_id")
    summ, allp = [], []
    for t in a.trees:
        clade = Path(t).parent.name
        tree = TreeNode.read(t)
        Dfull = tree.tip_tip_distances().to_data_frame()
        # Newick readers turn '_' in unquoted labels into ' '; sample ids never contain spaces (validated)
        Dfull.index = Dfull.columns = [str(i).replace(" ", "_") for i in Dfull.index]
        ids = [i for i in Dfull.index if i in meta.index]
        if len(ids) < 4:
            print(f"[strain] {clade}: <4 tips match the metadata; skipped", file=sys.stderr)
            continue
        D = Dfull.loc[ids, ids]
        iu = np.triu_indices(len(ids), 1)
        rows = pd.DataFrame({"clade": clade, "a": np.array(ids)[iu[0]], "b": np.array(ids)[iu[1]], "dist": D.to_numpy()[iu]})
        rows["subject_a"] = meta.loc[rows.a, "subject_id"].to_numpy()
        rows["subject_b"] = meta.loc[rows.b, "subject_id"].to_numpy()
        rows["within"] = rows.subject_a == rows.subject_b
        rows["diagnosis"] = np.where(rows.within, meta.loc[rows.a, "diagnosis"].to_numpy(), "between")
        rows["week_gap"] = np.where(rows.within, np.abs(meta.loc[rows.a, "week"].to_numpy() - meta.loc[rows.b, "week"].to_numpy()), np.nan)
        allp.append(rows)
        w, b = rows.loc[rows.within, "dist"], rows.loc[~rows.within, "dist"]
        if len(w) < 3:
            continue
        u = stats.mannwhitneyu(w, b, alternative="less")
        dx_groups = [g["dist"].to_numpy() for _, g in rows[rows.within].groupby("diagnosis") if len(g) >= 3]
        kw = stats.kruskal(*dx_groups).pvalue if len(dx_groups) > 1 else np.nan
        rho = stats.spearmanr(rows.loc[rows.within, "week_gap"], w) if w.nunique() > 1 else None
        summ.append({"clade": clade, "n_samples": len(ids), "n_subjects": meta.loc[ids, "subject_id"].nunique(),
                     "n_within_pairs": len(w), "median_within": w.median(), "median_between": b.median(),
                     "auc_within_lt_between": 1 - u.statistic / (len(w) * len(b)),
                     "p_within_lt_between": u.pvalue, "p_within_diff_by_diagnosis": kw,
                     "rho_within_dist_vs_week_gap": rho.statistic if rho is not None else np.nan})
    pd.DataFrame(summ).to_csv(out / "strain_summary.tsv", sep="\t", index=False)
    if allp:
        pd.concat(allp).to_csv(out / "strain_pairwise_distances.tsv.gz", sep="\t", index=False)
        d = pd.concat(allp)
        clades = d.clade.unique()
        fig, axes = plt.subplots(1, len(clades), figsize=(3.6 * len(clades) + 1, 4), squeeze=False)
        for ax, c in zip(axes[0], clades):
            sub = d[d.clade == c]
            order = ["between", "nonIBD", "CD", "UC"]
            data = [sub.loc[sub.diagnosis == k, "dist"] for k in order]
            bp = ax.boxplot(data, tick_labels=["between\nsubjects", "nonIBD", "CD", "UC"], patch_artist=True, showfliers=False)
            for patch, k in zip(bp["boxes"], order):
                patch.set_facecolor(COLORS.get(k, "#bbbbbb")); patch.set_alpha(0.8)
            ax.set_title(c, fontsize=10, fontweight="bold"); ax.set_ylabel("patristic distance")
        fig.suptitle("Strain persistence: same-subject pairs vs different-subject pairs", fontsize=10)
        fig.tight_layout()
        for ext in ("pdf", "png"):
            fig.savefig(Path(a.figdir) / f"Fig9_strain_persistence.{ext}", dpi=300)
    else:
        (out / "strain_summary.tsv").write_text("clade\tnote\nNA\tno StrainPhlAn clades configured\n")
    print(f"[strain] clades={len(summ)}", file=sys.stderr)


if __name__ == "__main__":
    main()

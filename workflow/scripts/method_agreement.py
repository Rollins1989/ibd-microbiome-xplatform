#!/usr/bin/env python3
"""Agreement between two shotgun profilers (MetaPhlAn vs Kraken2/Bracken) at species level."""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
import mbstats as mb  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metaphlan", required=True)
    ap.add_argument("--bracken", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    A, B = mb.read_layer(a.metaphlan), mb.read_layer(a.bracken)
    key = lambda idx: [mb.normalize_taxon_name(i.split(" ")[0]) + "_" + (i.split(" ", 1)[1].lower() if " " in i else "")  # noqa: E731
                       if mb.normalize_taxon_name(i.split(" ")[0]) else None for i in idx]
    A.index, B.index = key(A.index), key(B.index)
    A, B = A[A.index.notna()].groupby(level=0).sum(), B[B.index.notna()].groupby(level=0).sum()
    sp, sa = sorted(set(A.index) & set(B.index)), sorted(set(A.columns) & set(B.columns))
    if len(sp) < 5 or len(sa) < 5:
        raise SystemExit("too little overlap between profilers")
    Pa, Pb = mb.to_proportions(A.loc[sp, sa]), mb.to_proportions(B.loc[sp, sa])
    bc = (Pa - Pb).abs().sum(axis=0) / 2          # Bray-Curtis on closed shared-species profiles
    rows = []
    for s in sp:
        if Pa.loc[s].std() > 0 and Pb.loc[s].std() > 0:
            r = stats.spearmanr(Pa.loc[s], Pb.loc[s])
            rows.append({"species": s, "rho": r.statistic, "p": r.pvalue})
    sp_df = pd.DataFrame(rows)
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    sp_df.to_csv(out / "species_spearman_metaphlan_vs_bracken.tsv", sep="\t", index=False)
    pd.DataFrame([{"n_shared_species": len(sp), "n_samples": len(sa), "median_species_rho": sp_df.rho.median(),
                   "median_sample_bray_curtis": bc.median(),
                   "pct_metaphlan_abundance_shared": 100 * float(mb.to_proportions(A[sa]).loc[sp].sum().mean())}]
                 ).to_csv(out / "profiler_agreement_summary.tsv", sep="\t", index=False)
    print("[method_agreement] done", file=sys.stderr)


if __name__ == "__main__":
    main()

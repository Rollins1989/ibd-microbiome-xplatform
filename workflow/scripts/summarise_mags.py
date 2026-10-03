#!/usr/bin/env python3
"""Create the MAG catalogue (mag_summary.tsv) from dRep representatives + CheckM2 + GTDB-Tk.

Columns: mag, source_sample, completeness, contamination, quality, gtdb_classification,
classification_method, fastani_ani, genome_size
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rep-dir", required=True, help="dRep dereplicated_genomes directory")
    ap.add_argument("--quality", required=True, help="all_bins_quality.tsv from collect_bins.py")
    ap.add_argument("--gtdbtk-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    reps = sorted(p.name for p in Path(a.rep_dir).glob("*.fa"))
    q = pd.read_csv(a.quality, sep="\t").set_index("genome")
    frames = [pd.read_csv(f, sep="\t") for f in Path(a.gtdbtk_dir).rglob("gtdbtk.*.summary.tsv")]
    if not frames:
        raise SystemExit(f"no gtdbtk.*.summary.tsv under {a.gtdbtk_dir}")
    g = pd.concat(frames).set_index("user_genome")
    rows = []
    for fa in reps:
        stem = fa[:-3]
        gt = g.loc[stem] if stem in g.index else None
        rows.append({
            "mag": stem, "source_sample": q.loc[fa, "sample"], "completeness": q.loc[fa, "completeness"],
            "contamination": q.loc[fa, "contamination"], "quality": q.loc[fa, "quality"],
            "gtdb_classification": gt["classification"] if gt is not None else "Unclassified",
            "classification_method": gt.get("classification_method", np.nan) if gt is not None else np.nan,
            "fastani_ani": pd.to_numeric(gt.get("fastani_ani", np.nan), errors="coerce") if gt is not None else np.nan,
            "genome_size": q.loc[fa, "genome_size"]})
    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False, na_rep="NA")
    print(f"[summarise_mags] {len(rows)} representative MAGs", file=sys.stderr)


if __name__ == "__main__":
    main()

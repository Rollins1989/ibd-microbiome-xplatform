#!/usr/bin/env python3
"""Merge per-sample Bracken outputs into a species x sample table of re-estimated read counts."""
import argparse
import sys
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cols = {}
    for s, p in zip(a.samples, a.inputs):
        d = pd.read_csv(p, sep="\t")
        if "new_est_reads" not in d.columns or "name" not in d.columns:
            raise SystemExit(f"{p}: Bracken columns 'name'/'new_est_reads' not found: {list(d.columns)}")
        cols[s] = d.groupby("name")["new_est_reads"].sum()
    df = pd.DataFrame(cols).fillna(0.0)
    df.index.name = "feature"
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out, sep="\t")
    print(f"[merge_bracken] {df.shape[0]} species x {df.shape[1]} samples", file=sys.stderr)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Merge per-sample `coverm genome` tables into MAG x sample matrices (relative abundance %, covered fraction)."""
import argparse
import sys
from pathlib import Path

import pandas as pd


def pick(df, suffix, path):
    c = [x for x in df.columns if x.endswith(suffix)]
    if len(c) != 1:
        raise SystemExit(f"{path}: expected exactly one column ending in '{suffix}', found {c}")
    return df.set_index("Genome")[c[0]]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--out-abundance", required=True)
    ap.add_argument("--out-covered", required=True)
    a = ap.parse_args()
    ab, cv = {}, {}
    for s, p in zip(a.samples, a.inputs):
        d = pd.read_csv(p, sep="\t")
        ab[s] = pick(d, "Relative Abundance (%)", p)
        cv[s] = pick(d, "Covered Fraction", p)
    for df, path in ((pd.DataFrame(ab).fillna(0.0), a.out_abundance), (pd.DataFrame(cv).fillna(0.0), a.out_covered)):
        df.index.name = "mag"
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path, sep="\t")
    print(f"[merge_coverm] {len(a.samples)} samples", file=sys.stderr)


if __name__ == "__main__":
    main()

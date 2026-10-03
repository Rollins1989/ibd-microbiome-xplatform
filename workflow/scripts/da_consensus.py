#!/usr/bin/env python3
"""Combine per-method DA tables into a consensus call.

A feature/contrast is *consensus-significant* when q < threshold in >= min_methods methods
(capped at the number of methods that actually produced results) AND all significant methods
agree on the sign of the effect.
"""
import argparse
import sys

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--q", type=float, default=0.1)
    ap.add_argument("--min-methods", type=int, default=2)
    ap.add_argument("--primary", default="clr_lmm", help="method whose effect size is reported")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    df = pd.concat([pd.read_csv(p, sep="\t") for p in a.inputs], ignore_index=True)
    df["sig"] = df["qvalue"] < a.q
    df["sign"] = np.sign(df["estimate"])
    methods = sorted(df["method"].unique())
    need = min(a.min_methods, len(methods))
    rows = []
    for (layer, contrast, feat), g in df.groupby(["layer", "contrast", "feature"]):
        sg = g[g.sig]
        agree = sg["sign"].nunique() <= 1
        prim = g[g.method == a.primary]
        prim = prim if len(prim) else g.iloc[[0]]
        rows.append({"layer": layer, "contrast": contrast, "feature": feat, "n_methods": g["method"].nunique(),
                     "n_sig": len(sg), "methods_sig": ",".join(sorted(sg["method"])),
                     "direction": ("up" if sg["sign"].iloc[0] > 0 else "down") if len(sg) and agree else "mixed" if len(sg) else "ns",
                     "effect": prim["estimate"].iloc[0], "primary_method": prim["method"].iloc[0],
                     "min_q": g["qvalue"].min(),
                     "consensus": bool(len(sg) >= need and agree and need > 0)})
    out = pd.DataFrame(rows).sort_values(["layer", "contrast", "consensus", "n_sig", "min_q"],
                                         ascending=[True, True, False, False, True])
    out.to_csv(a.out, sep="\t", index=False)
    print(f"[da_consensus] methods={methods} need>={need} consensus features={int(out.consensus.sum())}", file=sys.stderr)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Merge per-sample MetaPhlAn 4 `-t rel_ab_w_read_stats` profiles.

Writes two clade x sample tables (all taxonomic levels, full lineage as row name):
  * relative abundance (%)
  * estimated number of reads from the clade
and a QC table with the UNCLASSIFIED fraction. Column positions are detected from the
header, and the script fails loudly if the expected columns are absent (format drift).
"""
import argparse
import sys
from pathlib import Path

import pandas as pd


def parse(path: Path) -> pd.DataFrame:
    header, rows = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                if line.startswith("#clade_name") or line.startswith("#taxon"):
                    header = line.lstrip("#").split("\t")
                continue
            rows.append(line.split("\t"))
    if header is None:
        raise SystemExit(f"{path}: no '#clade_name' header; was MetaPhlAn run with -t rel_ab_w_read_stats?")
    ra = [i for i, h in enumerate(header) if h == "relative_abundance"]
    rc = [i for i, h in enumerate(header) if h.startswith("estimated_number_of_reads")]
    if not ra or not rc:
        raise SystemExit(f"{path}: columns 'relative_abundance'/'estimated_number_of_reads...' not found in {header}")
    out = pd.DataFrame({"clade": [r[0] for r in rows],
                        "rel": [float(r[ra[0]]) for r in rows],
                        "reads": [float(r[rc[0]]) if len(r) > rc[0] else 0.0 for r in rows]})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--profiles", nargs="+", required=True)
    ap.add_argument("--samples", nargs="+", required=True, help="sample ids, same order as --profiles")
    ap.add_argument("--out-rel", required=True)
    ap.add_argument("--out-counts", required=True)
    ap.add_argument("--out-qc", required=True)
    a = ap.parse_args()
    if len(a.profiles) != len(a.samples):
        raise SystemExit("--profiles and --samples must have equal length")

    rel, cnt, qc = {}, {}, []
    for s, p in zip(a.samples, a.profiles):
        d = parse(Path(p))
        unc = d[d.clade == "UNCLASSIFIED"]
        qc.append({"sample_id": s, "unclassified_pct": float(unc.rel.sum()) if len(unc) else 0.0})
        d = d[d.clade.str.startswith(("k__", "d__"))]
        rel[s] = d.set_index("clade")["rel"]
        cnt[s] = d.set_index("clade")["reads"]
    for df, path in ((pd.DataFrame(rel).fillna(0.0), a.out_rel), (pd.DataFrame(cnt).fillna(0.0), a.out_counts)):
        df.index.name = "clade_name"
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path, sep="\t")
    pd.DataFrame(qc).to_csv(a.out_qc, sep="\t", index=False)
    print(f"[merge_metaphlan] {len(a.samples)} samples", file=sys.stderr)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Collect DAS Tool bins from all assembled samples, apply MIMAG quality thresholds and prepare dRep input.

Outputs
  <outdir>/passing/<sample>__<bin>.fa   bins with completeness >= min_comp and contamination <= max_cont
  <outdir>/genome_info.csv              genome,completeness,contamination   (dRep --genomeInfo format)
  <outdir>/all_bins_quality.tsv         every bin with CheckM2 metrics, sample and pass/fail
"""
import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--bin-dirs", nargs="+", required=True)
    ap.add_argument("--reports", nargs="+", required=True, help="CheckM2 quality_report.tsv per sample")
    ap.add_argument("--min-completeness", type=float, default=50)
    ap.add_argument("--max-contamination", type=float, default=10)
    ap.add_argument("--hq-completeness", type=float, default=90)
    ap.add_argument("--hq-contamination", type=float, default=5)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    if not (len(a.samples) == len(a.bin_dirs) == len(a.reports)):
        raise SystemExit("--samples, --bin-dirs and --reports must have the same length")
    out = Path(a.outdir)
    (out / "passing").mkdir(parents=True, exist_ok=True)
    rows = []
    for s, bd, rp in zip(a.samples, a.bin_dirs, a.reports):
        rp = Path(rp)
        if not rp.exists() or rp.stat().st_size == 0:
            continue
        q = pd.read_csv(rp, sep="\t")
        for _, r in q.iterrows():
            cands = [p for p in Path(bd).glob(f"{r['Name']}.*") if p.suffix in (".fa", ".fasta", ".fna")]
            if not cands:
                continue
            ok = r["Completeness"] >= a.min_completeness and r["Contamination"] <= a.max_contamination
            genome = f"{s}__{r['Name']}.fa"
            if ok:
                shutil.copy(cands[0], out / "passing" / genome)
            rows.append({"genome": genome, "sample": s, "bin": r["Name"], "completeness": r["Completeness"],
                         "contamination": r["Contamination"], "genome_size": r.get("Genome_Size"),
                         "quality": "HQ" if (r["Completeness"] >= a.hq_completeness and r["Contamination"] < a.hq_contamination)
                         else ("MQ" if ok else "low"), "pass": ok})
    df = pd.DataFrame(rows, columns=["genome", "sample", "bin", "completeness", "contamination", "genome_size", "quality", "pass"])
    df.to_csv(out / "all_bins_quality.tsv", sep="\t", index=False)
    df[df["pass"]][["genome", "completeness", "contamination"]].to_csv(out / "genome_info.csv", index=False)
    print(f"[collect_bins] {len(df)} bins, {int(df['pass'].sum())} pass (HQ={int((df.quality == 'HQ').sum())})", file=sys.stderr)


if __name__ == "__main__":
    main()

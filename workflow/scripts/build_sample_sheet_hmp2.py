#!/usr/bin/env python3
"""Build config/samples.tsv for HMP2/IBDMDB (BioProject PRJNA398089) from two public tables.

Inputs
  --metadata  HMP2 sample metadata CSV/TSV from https://ibdmdb.org (see docs/databases.md)
  --ena       ENA read-run report, fetched with
              curl -o ena.tsv 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA398089&result=read_run&fields=run_accession,sample_alias,experiment_alias,library_strategy,library_source,instrument_model,first_created,fastq_ftp,fastq_md5&format=tsv'
Behaviour
  * joins on the configured specimen identifier, reports the join rate and every dropped record (reason + count)
  * keeps only paired-end runs; one sample per (specimen, platform)
  * `batch` is a PROXY (ENA submission month + instrument) unless --batch-column names a metadata column
  * never invents values: unmapped diagnoses/platforms are dropped and counted
"""
import argparse
import sys
from pathlib import Path

import pandas as pd
import yaml


def read_any(path):
    sep = "\t" if str(path).endswith((".tsv", ".txt")) else ","
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, na_values=["", "NA", "N/A", "nan"])


def need(df, col, what):
    if col not in df.columns:
        raise SystemExit(f"column '{col}' ({what}) not found. Available: {sorted(df.columns)}\n"
                         "-> adjust config/hmp2_columns.yaml")
    return df[col]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--ena", required=True)
    ap.add_argument("--mapping", default="config/hmp2_columns.yaml")
    ap.add_argument("--filter", action="append", default=[], metavar="COL=V1,V2",
                    help="keep rows whose metadata COL is in the list (repeatable), e.g. biopsy/stool column")
    ap.add_argument("--batch-column", default=None)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.mapping))
    mm = cfg["metadata"]

    meta = read_any(a.metadata)
    ena = read_any(a.ena)
    log = []
    for f in a.filter:
        col, vals = f.split("=", 1)
        n0 = len(meta)
        meta = meta[need(meta, col, "filter").isin(vals.split(","))]
        log.append((f"filter {f}", n0 - len(meta)))

    key = need(meta, mm["join_key"], "specimen id")
    out = pd.DataFrame({"specimen_id": key.astype(str)})
    out["subject_id"] = need(meta, mm["subject"], "subject").astype(str)
    out["data_type"] = need(meta, mm["data_type"], "data type").map(cfg["data_type_map"])
    out["week"] = pd.to_numeric(need(meta, mm["week"], "week"), errors="coerce")
    out["diagnosis"] = need(meta, mm["diagnosis"], "diagnosis").map(cfg["diagnosis_map"])
    for tgt in ("age", "sex", "bmi", "antibiotics", "immunosuppressants", "hbi", "sccai"):
        src = mm.get(tgt)
        out[tgt] = meta[src].to_numpy() if src else pd.NA
        if src and src not in meta.columns:
            raise SystemExit(f"column '{src}' ({tgt}) not found. Available: {sorted(meta.columns)}")
    out["_batch_meta"] = meta[a.batch_column].to_numpy() if a.batch_column else pd.NA
    for col, why in (("data_type", "unmapped platform"), ("diagnosis", "unmapped diagnosis"), ("week", "missing week")):
        bad = out[col].isna()
        log.append((why, int(bad.sum())))
        out = out[~bad]

    ena = ena.rename(columns={cfg["metadata"]["ena_join_field"]: "_key"})
    need(ena, "_key", "ENA join field")
    ena = ena.dropna(subset=["fastq_ftp"])
    ena = ena[ena.fastq_ftp.str.count(";") == 1]               # exactly two files = paired-end
    merged = out.merge(ena, left_on="specimen_id", right_on="_key", how="left", indicator=True)
    log.append(("metadata records without an ENA run", int((merged["_merge"] == "left_only").sum())))
    merged = merged[merged["_merge"] == "both"].copy()
    # platform consistency with ENA library strategy
    exp = merged.data_type.map(cfg["ena_library_strategy"])
    mism = exp != merged.library_strategy
    log.append(("platform/library_strategy mismatch (dropped)", int(mism.sum())))
    merged = merged[~mism]
    merged = merged.sort_values("run_accession").drop_duplicates(["specimen_id", "data_type"], keep="first")

    urls = merged.fastq_ftp.str.split(";", expand=True)
    sid = merged.specimen_id + "_" + merged.data_type
    batch = merged["_batch_meta"].where(merged["_batch_meta"].notna(),
                                        merged.first_created.str[:7].fillna("NA") + "|" + merged.instrument_model.fillna("NA"))
    batch = batch.str.replace(r"[^A-Za-z0-9_.\-]+", "_", regex=True)
    res = pd.DataFrame({
        "sample_id": sid.str.replace(r"[^A-Za-z0-9_.\-]+", "_", regex=True), "specimen_id": merged.specimen_id,
        "subject_id": merged.subject_id, "data_type": merged.data_type, "batch": batch, "sample_type": "sample",
        "week": merged.week, "diagnosis": merged.diagnosis, "age": merged.age, "sex": merged.sex, "bmi": merged.bmi,
        "antibiotics": merged.antibiotics, "immunosuppressants": merged.immunosuppressants, "hbi": merged.hbi,
        "sccai": merged.sccai, "assemble": "no", "fastq_1": "https://" + urls[0], "fastq_2": "https://" + urls[1]})
    # one diagnosis per subject (HMP2 occasionally reassigns): keep the most frequent
    dx = res.groupby("subject_id")["diagnosis"].agg(lambda s: s.value_counts().index[0])
    changed = int((res.diagnosis != res.subject_id.map(dx)).sum())
    res["diagnosis"] = res.subject_id.map(dx)
    log.append(("records whose diagnosis was harmonised per subject", changed))
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(a.out, sep="\t", index=False, na_rep="NA")
    print("== sample sheet build report ==", file=sys.stderr)
    for k, v in log:
        print(f"  {k}: {v}", file=sys.stderr)
    print(f"  written: {len(res)} samples, {res.subject_id.nunique()} subjects, "
          f"{res.groupby('specimen_id').data_type.nunique().eq(2).sum()} paired specimens -> {a.out}", file=sys.stderr)
    print("  NOTE: set `assemble` to yes for the subset of shotgun samples to assemble; review `batch` (proxy).", file=sys.stderr)


if __name__ == "__main__":
    main()

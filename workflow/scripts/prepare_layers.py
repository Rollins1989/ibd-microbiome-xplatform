#!/usr/bin/env python3
"""Build the standardised feature tables ("layers") used by every downstream analysis.

Every layer is a features x samples TSV (first column `feature`). Samples failing the depth
filters are removed from *all* layers of their platform and reported in layer_summary.tsv.
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
import mbstats as mb  # noqa: E402

RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]


def genus_labels(tax: pd.DataFrame) -> pd.Series:
    """Genus label; unresolved genera become 'Unclassified_<lowest assigned rank name>'."""
    out = {}
    for asv, r in tax.iterrows():
        if pd.notna(r.get("Genus")):
            out[asv] = str(r["Genus"])
            continue
        low = next((str(r[k]) for k in reversed(RANKS[:-1]) if k in r and pd.notna(r[k])), "Bacteria")
        out[asv] = f"Unclassified_{low}"
    return pd.Series(out)


def mpa_rank(df: pd.DataFrame, prefix: str, integer: bool) -> pd.DataFrame:
    last = df.index.to_series().str.split("|").str[-1]
    sel = last.str.startswith(prefix) & ~df.index.to_series().str.contains(r"\|t__")
    sub = df.loc[sel.values].copy()
    sub.index = [i.split("|")[-1][len(prefix):].replace("_", " ") for i in sub.index]
    sub = sub.groupby(level=0).sum()
    return sub.round().astype(float) if integer else sub


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--samples-meta", required=True)
    ap.add_argument("--specimens", required=True)
    ap.add_argument("--asv", required=True)
    ap.add_argument("--taxonomy", required=True)
    ap.add_argument("--mpa-counts", required=True)
    ap.add_argument("--humann", default=None)
    ap.add_argument("--picrust", required=True)
    ap.add_argument("--mag-abundance", default=None)
    ap.add_argument("--min-reads-16s", type=float, default=5000)
    ap.add_argument("--min-reads-mgx", type=float, default=100000)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(a.samples_meta, sep="\t", na_values=mb.NA_STRINGS, keep_default_na=False)
    spec = pd.read_csv(a.specimens, sep="\t", na_values=mb.NA_STRINGS, keep_default_na=False)
    summary, dropped = [], []

    # ------------------------------------------------------------- 16S
    asv = mb.read_layer(a.asv)
    tax = pd.read_csv(a.taxonomy, sep="\t", index_col=0, na_values=mb.NA_STRINGS, keep_default_na=False)
    ids16 = [s for s in meta.loc[meta.data_type == "16S", "sample_id"] if s in asv.columns]
    asv = asv[ids16]
    depth = asv.sum(axis=0)
    bad = depth[depth < a.min_reads_16s].index.tolist()
    dropped += [(s, "16S", f"depth<{a.min_reads_16s:g}") for s in bad]
    asv = asv.drop(columns=bad)
    asv = asv.loc[asv.sum(axis=1) > 0]
    mb.write_layer(asv, out / "asv_table.filtered.tsv", "ASV")
    lab = genus_labels(tax.loc[asv.index.intersection(tax.index)])
    g16 = asv.loc[lab.index].groupby(lab.values).sum()
    mb.write_layer(g16, out / "16s_genus.tsv")
    pic = mb.read_layer(a.picrust)
    pic = pic[[c for c in pic.columns if c in g16.columns]]
    mb.write_layer(pic, out / "16s_pathway.tsv")
    keep16 = set(g16.columns)

    # ------------------------------------------------------------- MGX
    mpa = mb.read_layer(a.mpa_counts)
    idsm = [s for s in meta.loc[meta.data_type == "MGX", "sample_id"] if s in mpa.columns]
    mpa = mpa[idsm]
    sp = mpa_rank(mpa, "s__", True)
    depth = sp.sum(axis=0)
    bad = depth[depth < a.min_reads_mgx].index.tolist()
    dropped += [(s, "MGX", f"microbial reads<{a.min_reads_mgx:g}") for s in bad]
    mpa = mpa.drop(columns=bad)
    sp = mpa_rank(mpa, "s__", True)
    sp = sp.loc[sp.sum(axis=1) > 0]
    ge = mpa_rank(mpa, "g__", True)
    ge = ge.loc[ge.sum(axis=1) > 0]
    mb.write_layer(sp, out / "mgx_species.tsv")
    mb.write_layer(ge, out / "mgx_genus.tsv")
    keepm = set(sp.columns)
    if a.humann and Path(a.humann).exists():
        h = pd.read_csv(a.humann, sep="\t", index_col=0)
        h.index = h.index.astype(str)
        h.columns = [re.sub(r"_Abundance(-CPM|-RPKs)?$", "", c) for c in h.columns]
        h = h.loc[~h.index.isin(["UNMAPPED", "UNINTEGRATED"])]
        h = h.loc[:, [c for c in h.columns if c in keepm]]
        mb.write_layer(h.astype(float), out / "mgx_pathway.tsv")
    if a.mag_abundance and Path(a.mag_abundance).exists():
        m = mb.read_layer(a.mag_abundance)
        m = m.loc[~m.index.isin(["unmapped", "Unmapped"])]
        m = m.loc[:, [c for c in m.columns if c in keepm]]
        mb.write_layer(m, out / "mag.tsv")

    # ------------------------------------------------------------- metadata
    keep = keep16 | keepm
    mfin = meta[meta.sample_id.isin(keep)].copy()
    mfin.to_csv(out / "metadata.analysis.tsv", sep="\t", index=False, na_rep="NA")
    for col, plat_keep in (("sample_16s", keep16), ("sample_mgx", keepm)):
        if col in spec:
            spec[col] = spec[col].where(spec[col].isin(plat_keep))
    spec["paired"] = spec[["sample_16s", "sample_mgx"]].notna().all(axis=1)
    spec = spec[spec[["sample_16s", "sample_mgx"]].notna().any(axis=1)]
    spec.to_csv(out / "specimens.analysis.tsv", sep="\t", index=False, na_rep="NA")

    for f in sorted(out.glob("*.tsv")):
        if f.name.startswith(("16s_", "mgx_", "mag")) and "asv" not in f.name:
            d = mb.read_layer(f)
            summary.append({"layer": f.stem, "features": d.shape[0], "samples": d.shape[1]})
    pd.DataFrame(summary).to_csv(out / "layer_summary.tsv", sep="\t", index=False)
    pd.DataFrame(dropped, columns=["sample_id", "platform", "reason"]).to_csv(out / "dropped_samples.tsv", sep="\t", index=False)
    print(f"[prepare_layers] retained 16S={len(keep16)} MGX={len(keepm)}; dropped={len(dropped)}", file=sys.stderr)


if __name__ == "__main__":
    main()

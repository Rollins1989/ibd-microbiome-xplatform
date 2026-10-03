#!/usr/bin/env python3
"""Simulate upstream outputs (formats identical to the real workflow) with planted signal.

Planted truth (so tests can verify recovery):
  * Escherichia up, Faecalibacterium & Roseburia down in CD; Escherichia + Ruminococcus up in UC
  * activity: Escherichia up, Faecalibacterium down
  * pre-flare: Veillonella up (weak)
  * 16S and shotgun share one latent community per specimen (+ platform bias/noise)
  * 3 contaminant ASVs present in negative controls
Usage: python tests/simulate_data.py tests/sim [--seed 7]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

TAX = [  # genus, family, order, class, phylum, species list
    ("Faecalibacterium", "Oscillospiraceae", "Eubacteriales", "Clostridia", "Firmicutes", ["prausnitzii", "sp_A"]),
    ("Roseburia", "Lachnospiraceae", "Eubacteriales", "Clostridia", "Firmicutes", ["intestinalis", "hominis"]),
    ("Blautia", "Lachnospiraceae", "Eubacteriales", "Clostridia", "Firmicutes", ["obeum", "wexlerae"]),
    ("Ruminococcus", "Ruminococcaceae", "Eubacteriales", "Clostridia", "Firmicutes", ["gnavus", "bromii"]),
    ("Eubacterium", "Eubacteriaceae", "Eubacteriales", "Clostridia", "Firmicutes", ["rectale", "hallii"]),
    ("Anaerostipes", "Lachnospiraceae", "Eubacteriales", "Clostridia", "Firmicutes", ["hadrus"]),
    ("Coprococcus", "Lachnospiraceae", "Eubacteriales", "Clostridia", "Firmicutes", ["comes", "catus"]),
    ("Dorea", "Lachnospiraceae", "Eubacteriales", "Clostridia", "Firmicutes", ["longicatena"]),
    ("Veillonella", "Veillonellaceae", "Veillonellales", "Negativicutes", "Firmicutes", ["parvula", "atypica"]),
    ("Streptococcus", "Streptococcaceae", "Lactobacillales", "Bacilli", "Firmicutes", ["salivarius", "parasanguinis"]),
    ("Lactobacillus", "Lactobacillaceae", "Lactobacillales", "Bacilli", "Firmicutes", ["rhamnosus"]),
    ("Enterococcus", "Enterococcaceae", "Lactobacillales", "Bacilli", "Firmicutes", ["faecalis"]),
    ("Bacteroides", "Bacteroidaceae", "Bacteroidales", "Bacteroidia", "Bacteroidota", ["vulgatus", "fragilis", "uniformis"]),
    ("Prevotella", "Prevotellaceae", "Bacteroidales", "Bacteroidia", "Bacteroidota", ["copri"]),
    ("Alistipes", "Rikenellaceae", "Bacteroidales", "Bacteroidia", "Bacteroidota", ["putredinis", "onderdonkii"]),
    ("Parabacteroides", "Tannerellaceae", "Bacteroidales", "Bacteroidia", "Bacteroidota", ["distasonis"]),
    ("Bifidobacterium", "Bifidobacteriaceae", "Bifidobacteriales", "Actinobacteria", "Actinobacteria", ["longum", "adolescentis"]),
    ("Collinsella", "Coriobacteriaceae", "Coriobacteriales", "Coriobacteriia", "Actinobacteria", ["aerofaciens"]),
    ("Escherichia", "Enterobacteriaceae", "Enterobacterales", "Gammaproteobacteria", "Proteobacteria", ["coli"]),
    ("Klebsiella", "Enterobacteriaceae", "Enterobacterales", "Gammaproteobacteria", "Proteobacteria", ["pneumoniae"]),
    ("Akkermansia", "Akkermansiaceae", "Verrucomicrobiales", "Verrucomicrobiae", "Verrucomicrobiota", ["muciniphila"]),
    ("Fusobacterium", "Fusobacteriaceae", "Fusobacteriales", "Fusobacteriia", "Fusobacteriota", ["nucleatum"]),
    ("Haemophilus", "Pasteurellaceae", "Pasteurellales", "Gammaproteobacteria", "Proteobacteria", ["parainfluenzae"]),
    ("Sutterella", "Sutterellaceae", "Burkholderiales", "Betaproteobacteria", "Proteobacteria", ["wadsworthensis"]),
    ("Desulfovibrio", "Desulfovibrionaceae", "Desulfovibrionales", "Desulfovibrionia", "Desulfobacterota", ["piger"]),
    ("Clostridium", "Clostridiaceae", "Eubacteriales", "Clostridia", "Firmicutes", ["innocuum", "symbiosum"]),
    ("Oscillibacter", "Oscillospiraceae", "Eubacteriales", "Clostridia", "Firmicutes", ["sp_B"]),
    ("Subdoligranulum", "Oscillospiraceae", "Eubacteriales", "Clostridia", "Firmicutes", ["variabile"]),
    ("Dialister", "Veillonellaceae", "Veillonellales", "Negativicutes", "Firmicutes", ["invisus"]),
    ("Methanobrevibacter", "Methanobacteriaceae", "Methanobacteriales", "Methanobacteria", "Euryarchaeota", ["smithii"]),
]
SILVA = {"Escherichia": "Escherichia-Shigella", "Eubacterium": "[Eubacterium] rectale group",
         "Ruminococcus": "[Ruminococcus] gnavus group", "Clostridium": "Clostridium sensu stricto 1"}
G = [t[0] for t in TAX]
NG = len(G)
AMR = ["acrF", "mdtM", "emrD", "tetQ", "ermG", "CfxA6", "aph(3')-IIa", "blaTEM-1", "vanXY", "msrC"]


def newick(node):
    """node = (label:str, branch_length) for a leaf, or (list_of_nodes, branch_length) for an internal node."""
    first, bl = node
    if isinstance(first, str):
        return f"{first}:{bl:.5f}"
    return "(" + ",".join(newick(k) for k in first) + f"):{bl:.5f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    out = Path(a.outdir)
    res = out / "results"
    for d in ["16s/dada2", "16s/tree", "16s/picrust2", "shotgun/metaphlan", "shotgun/humann",
              "shotgun/strainphlan/t__SGB0001", "mags/summary", "mags/tree", "mags/abundance"]:
        (res / d).mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ design
    subj = ([("N%02d" % i, "nonIBD") for i in range(20)] + [("C%02d" % i, "CD") for i in range(22)]
            + [("U%02d" % i, "UC") for i in range(22)])
    gi = {g: i for i, g in enumerate(G)}
    mu0 = rng.normal(0, 1.4, NG)
    dx_eff = {"nonIBD": np.zeros(NG), "CD": np.zeros(NG), "UC": np.zeros(NG)}
    for g, v in {"Escherichia": 2.2, "Fusobacterium": 1.5, "Faecalibacterium": -1.8, "Roseburia": -1.4, "Klebsiella": 1.0}.items():
        dx_eff["CD"][gi[g]] = v
    for g, v in {"Escherichia": 1.3, "Ruminococcus": 1.4, "Faecalibacterium": -1.1, "Akkermansia": -1.0}.items():
        dx_eff["UC"][gi[g]] = v
    act_eff = np.zeros(NG); act_eff[gi["Escherichia"]] = 0.9; act_eff[gi["Faecalibacterium"]] = -0.6
    pre_eff = np.zeros(NG); pre_eff[gi["Veillonella"]] = 0.9
    bias16 = rng.normal(0, 0.5, NG)
    # planted modules: Bacteroides/Alistipes/Parabacteroides co-vary (+); Bifidobacterium/Collinsella co-vary (+)
    # and the second module is anti-correlated with Prevotella
    MOD1 = np.zeros(NG); MOD2 = np.zeros(NG)
    for g in ("Bacteroides", "Alistipes", "Parabacteroides"):
        MOD1[gi[g]] = 1.0
    MOD1[gi["Prevotella"]] = -0.9
    for g in ("Bifidobacterium", "Collinsella", "Lactobacillus"):
        MOD2[gi[g]] = 1.0

    rows, latent = [], {}
    for sid, dx in subj:
        n_vis = rng.integers(5, 8)
        weeks = np.cumsum(np.r_[0, rng.choice([2, 3, 4], n_vis - 1)])
        age = float(rng.integers(18, 70)); sex = rng.choice(["F", "M"]); bmi = round(float(rng.normal(25, 4)), 1)
        u = rng.normal(0, 1.0, NG)  # subject-specific offset
        state = 0
        acts = []
        for w in weeks:
            if dx != "nonIBD":
                state = int(rng.random() < (0.7 if state else 0.3))
            acts.append(state)
        for k, (w, act) in enumerate(zip(weeks, acts)):
            nxt = acts[k + 1] if k + 1 < len(acts) else None
            pre = (act == 0 and nxt == 1)
            spec = f"{sid}_v{k}"
            f1, f2 = rng.normal(0, 1.3, 2)   # latent co-abundance factors (planted network modules)
            abu = (mu0 + dx_eff[dx] + u + act * act_eff + pre * pre_eff + rng.normal(0, 0.5, NG)
                   + f1 * MOD1 + f2 * MOD2)
            latent[spec] = abu
            hbi = sccai = np.nan
            if dx == "CD":
                hbi = int(rng.integers(5, 12) if act else rng.integers(0, 4))
            if dx == "UC":
                sccai = int(rng.integers(5, 12) if act else rng.integers(0, 4))
            rows.append(dict(specimen_id=spec, subject_id=sid, week=int(w), diagnosis=dx, age=age, sex=sex, bmi=bmi,
                             antibiotics=rng.choice(["Yes", "No"], p=[0.12, 0.88]),
                             immunosuppressants=rng.choice(["Yes", "No"], p=[0.3, 0.7]) if dx != "nonIBD" else "No",
                             hbi=hbi, sccai=sccai))
    spec_df = pd.DataFrame(rows)

    sample_rows = []
    for r in rows:
        for dt, p in (("16S", 0.85), ("MGX", 0.85)):
            if rng.random() < p:
                d = dict(r)
                d.update(sample_id=f"{r['specimen_id']}_{dt}", data_type=dt, sample_type="sample",
                         batch=rng.choice(["r1", "r2", "r3"]) if dt == "16S" else rng.choice(["b1", "b2"]),
                         assemble="yes" if dt == "MGX" and rng.random() < 0.3 else "no",
                         fastq_1="x_1.fastq.gz", fastq_2="x_2.fastq.gz")
                sample_rows.append(d)
    for i in range(6):
        sample_rows.append(dict(sample_id=f"NEGCTRL{i+1}", specimen_id=f"NEGCTRL{i+1}", subject_id=np.nan, data_type="16S",
                                sample_type="negative_control", batch="r1", week=np.nan, diagnosis=np.nan,
                                assemble="no", fastq_1="x_1.fastq.gz", fastq_2="x_2.fastq.gz"))
    sdf = pd.DataFrame(sample_rows)
    cols = ["sample_id", "specimen_id", "subject_id", "data_type", "batch", "sample_type", "week", "diagnosis", "age",
            "sex", "bmi", "antibiotics", "immunosuppressants", "hbi", "sccai", "assemble", "fastq_1", "fastq_2"]
    sdf = sdf[cols]
    sdf.to_csv(out / "samples.tsv", sep="\t", index=False, na_rep="NA")

    def softmax(x):
        e = np.exp(x - x.max()); return e / e.sum()

    # ----------------------------------------------------------------- 16S
    asv_meta = []
    for g_i, (g, fam, *_r) in enumerate(TAX):
        for k in range(rng.integers(2, 6)):
            asv_meta.append((f"ASV{len(asv_meta)+1:04d}", g_i, rng.dirichlet(np.ones(1) * 1)[0]))
    n_asv = len(asv_meta)
    asv_ids = [a[0] for a in asv_meta]
    asv_genus = np.array([a[1] for a in asv_meta])
    split = np.zeros(n_asv)
    for g_i in range(NG):
        idx = np.where(asv_genus == g_i)[0]
        split[idx] = rng.dirichlet(np.ones(len(idx)) * 0.8)
    contam = ["ASV9001", "ASV9002", "ASV9003"]
    s16 = sdf[sdf.data_type == "16S"]
    tab16 = pd.DataFrame(0, index=asv_ids + contam, columns=s16.sample_id)
    for _, r in s16.iterrows():
        if r.sample_type == "negative_control":
            tab16.loc[contam, r.sample_id] = rng.poisson(400, 3)
            tab16.loc[asv_ids, r.sample_id] = rng.poisson(0.3, n_asv)
            continue
        p = softmax(latent[r.specimen_id] + bias16)
        pa = p[asv_genus] * split
        depth = int(rng.lognormal(np.log(25000), 0.4))
        tab16.loc[asv_ids, r.sample_id] = rng.multinomial(depth, pa / pa.sum())
        if rng.random() < 0.3:
            tab16.loc[contam, r.sample_id] = rng.poisson(15, 3)
    tab16.index.name = "ASV"
    tab16.to_csv(res / "16s/dada2/asv_table.tsv", sep="\t")

    tax_rows = []
    for aid, g_i, _ in asv_meta:
        g, fam, order, cls, phy, sp = TAX[g_i]
        genus = None if rng.random() < 0.06 else SILVA.get(g, g)
        tax_rows.append([aid, "Bacteria" if phy != "Euryarchaeota" else "Archaea", phy, cls, order, fam, genus, None])
    for c in contam:
        tax_rows.append([c, "Bacteria", "Proteobacteria", "Gammaproteobacteria", "Burkholderiales", "Burkholderiaceae", "Ralstonia", None])
    pd.DataFrame(tax_rows, columns=["ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]).to_csv(
        res / "16s/dada2/taxonomy.tsv", sep="\t", index=False, na_rep="NA")

    kids = []
    for g_i in range(NG):
        leaves = [(asv_ids[i], rng.uniform(0.005, 0.03)) for i in np.where(asv_genus == g_i)[0]]
        kids.append((leaves, rng.uniform(0.05, 0.2)))
    kids.append(([(c, rng.uniform(0.05, 0.1)) for c in contam], 0.3))
    (res / "16s/tree/asv_tree.rooted.nwk").write_text(newick((kids, 0.0)) + ";\n")

    # ----------------------------------------------------------------- MGX
    sp_names, sp_genus, sp_split = [], [], []
    for g_i, t in enumerate(TAX):
        fr = rng.dirichlet(np.ones(len(t[5])) * 1.5)
        for s, f in zip(t[5], fr):
            sp_names.append(f"{t[0]}_{s}"); sp_genus.append(g_i); sp_split.append(f)
    sp_genus = np.array(sp_genus); sp_split = np.array(sp_split)
    smg = sdf[sdf.data_type == "MGX"]
    rel_sp = pd.DataFrame(0.0, index=sp_names, columns=smg.sample_id)
    cnt_sp = rel_sp.copy()
    for _, r in smg.iterrows():
        p = softmax(latent[r.specimen_id] + rng.normal(0, 0.25, NG))
        ps = p[sp_genus] * sp_split
        ps /= ps.sum()
        depth = int(rng.lognormal(np.log(4e6), 0.4))
        c = rng.multinomial(depth, ps)
        cnt_sp[r.sample_id] = c
        rel_sp[r.sample_id] = 100 * c / c.sum()

    def lineage(g_i, level):
        g, fam, order, cls, phy, _sp = TAX[g_i]
        k = "k__Archaea" if phy == "Euryarchaeota" else "k__Bacteria"
        parts = [k, f"p__{phy}", f"c__{cls}", f"o__{order}", f"f__{fam}", f"g__{g}"]
        return "|".join(parts[:level])

    def hier(df):
        recs = {}
        for name, g_i in zip(df.index, sp_genus):
            sp = name.split("_", 1)[1] if False else name
            for lvl in range(1, 7):
                recs.setdefault(lineage(g_i, lvl), 0)
                recs[lineage(g_i, lvl)] = recs[lineage(g_i, lvl)] + df.loc[name].to_numpy()
            recs[lineage(g_i, 6) + f"|s__{name}"] = df.loc[name].to_numpy()
        out_df = pd.DataFrame(recs, index=df.columns).T
        out_df.index.name = "clade_name"
        return out_df
    hier(rel_sp).to_csv(res / "shotgun/metaphlan/merged_relabund.tsv", sep="\t")
    hier(cnt_sp).to_csv(res / "shotgun/metaphlan/merged_estcounts.tsv", sep="\t")

    # pathways: sparse loadings on genus proportions
    npw = 40
    L16 = np.zeros((npw, NG))
    for j in range(npw):
        L16[j, rng.choice(NG, rng.integers(2, 6), replace=False)] = rng.uniform(0.5, 2, )
    pw = [f"PWY-{1000+j}: simulated pathway {j}" for j in range(npw)]
    for tag, sm, noise, fname, hdr in (("mgx", smg, 0.15, "shotgun/humann/merged_pathabundance_cpm_unstratified.tsv", "# Pathway"),
                                       ("16s", s16[s16.sample_type == "sample"], 0.35, "16s/picrust2/path_abun_unstrat.tsv", "pathway")):
        mat = pd.DataFrame(0.0, index=pw, columns=sm.sample_id)
        for _, r in sm.iterrows():
            p = softmax(latent[r.specimen_id] + (bias16 if tag == "16s" else 0))
            v = L16 @ p * np.exp(rng.normal(0, noise, npw))
            mat[r.sample_id] = v / v.sum() * (1e6 if tag == "mgx" else 5e5)
        if tag == "mgx":
            mat.loc["UNMAPPED"] = 3e5; mat.loc["UNINTEGRATED"] = 2e5
            mat.columns = [c + "_Abundance-CPM" for c in mat.columns]
        mat.index.name = hdr
        mat.to_csv(res / fname, sep="\t")

    # ----------------------------------------------------------------- MAGs
    n_mag = 36
    mag_sp = rng.choice(len(sp_names), n_mag, replace=False)
    ecoli = sp_names.index("Escherichia_coli")
    if ecoli not in mag_sp:
        mag_sp[0] = ecoli
    mags, tips = [], []
    for m, s in enumerate(mag_sp):
        g_i = sp_genus[s]; g, fam, order, cls, phy, _ = TAX[g_i]
        novel = (rng.random() < 0.25) or (m == 1)
        sp_label = "s__" if novel else f"s__{sp_names[s].replace('_', ' ')}"
        name = f"MAG{m+1:03d}"
        comp = rng.uniform(55, 99); con = rng.uniform(0.2, 9)
        mags.append(dict(mag=name, source_sample=rng.choice(smg.sample_id), completeness=comp, contamination=con,
                         quality="HQ" if comp >= 90 and con < 5 else "MQ",
                         gtdb_classification=f"d__Bacteria;p__{phy};c__{cls};o__{order};f__{fam};g__{g};{sp_label}",
                         classification_method="taxonomic placement" if novel else "ANI",
                         fastani_ani=np.nan if novel else rng.uniform(95.2, 99.5), genome_size=int(rng.uniform(1.8e6, 5.5e6)),
                         ref_species=sp_names[s]))
        tips.append(name)
    mags = pd.DataFrame(mags)
    mags.to_csv(res / "mags/summary/mag_summary.tsv", sep="\t", index=False, na_rep="NA")
    ab = pd.DataFrame(0.0, index=mags.mag, columns=smg.sample_id)
    for _, r in smg.iterrows():
        p = softmax(latent[r.specimen_id] + rng.normal(0, 0.25, NG))
        v = np.array([p[sp_genus[s]] * sp_split[s] for s in mag_sp])
        v = v / v.sum() * rng.uniform(40, 75)
        v[v < 0.02] = 0
        ab[r.sample_id] = v
    ab.loc["unmapped"] = 100 - ab.sum()
    ab.index.name = "mag"
    ab.to_csv(res / "mags/abundance/mag_relabund.tsv", sep="\t")
    amr = []
    for _, m in mags.iterrows():
        g = m.ref_species.split("_")[0]
        n = rng.integers(2, 7) if g in ("Escherichia", "Klebsiella", "Enterococcus") else rng.integers(0, 2)
        for gene in rng.choice(AMR, n, replace=False):
            amr.append(dict(mag=m.mag, source="rgi", gene=gene, drug_class="simulated", identity=round(rng.uniform(80, 100), 1)))
    pd.DataFrame(amr, columns=["mag", "source", "gene", "drug_class", "identity"]).to_csv(
        res / "mags/summary/amr_hits.tsv", sep="\t", index=False)
    grp = {}
    for t, s in zip(tips, mag_sp):
        grp.setdefault(sp_genus[s], []).append((t, rng.uniform(0.02, 0.1)))
    (res / "mags/tree/mag_tree.nwk").write_text(
        newick(([(v, rng.uniform(0.1, 0.3)) for v in grp.values()], 0.0)) + ";\n")

    # ------------------------------------------------- StrainPhlAn-like tree
    by_subj = smg.groupby("subject_id")["sample_id"].apply(list)
    clades = []
    for s, ids in by_subj.items():
        clades.append(([(i, rng.uniform(0.0005, 0.004)) for i in ids], rng.uniform(0.05, 0.25)))
    (res / "shotgun/strainphlan/t__SGB0001/tree.nwk").write_text(newick((clades, 0.0)) + ";\n")
    print(f"simulated {len(spec_df)} specimens, {len(sdf)} samples -> {out}")


if __name__ == "__main__":
    main()

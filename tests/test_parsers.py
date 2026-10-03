"""Glue-script tests on tiny files that mimic the real tool output formats (CheckM2, GTDB-Tk, RGI, ABRicate, CoverM,
MetaPhlAn rel_ab_w_read_stats, Bracken). They guard the parsing/merging code; they cannot guard against a
tool changing its format - the parsers raise explicit errors in that case."""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
S = ROOT / "workflow" / "scripts"


def run(script, *args):
    r = subprocess.run([sys.executable, str(S / script), *map(str, args)], capture_output=True, text=True)
    return r


def test_merge_metaphlan(tmp_path):
    hdr = "#mpa_vJan21_CHOCOPhlAnSGB_202103\n#clade_name\tNCBI_tax_id\trelative_abundance\tcoverage\testimated_number_of_reads_from_the_clade\n"
    p1 = hdr + ("UNCLASSIFIED\t-1\t20.0\t0\t0\n" "k__Bacteria\t2\t80.0\t5\t800\n"
                "k__Bacteria|p__Firmicutes|c__C|o__O|f__F|g__Blautia\t1|2|3|4|5|6\t80.0\t5\t800\n"
                "k__Bacteria|p__Firmicutes|c__C|o__O|f__F|g__Blautia|s__Blautia_obeum\t1|2|3|4|5|6|7\t80.0\t5\t800\n")
    (tmp_path / "a.tsv").write_text(p1)
    (tmp_path / "b.tsv").write_text(p1.replace("800", "400").replace("80.0", "100.0").replace("20.0", "0.0"))
    r = run("merge_metaphlan.py", "--profiles", tmp_path / "a.tsv", tmp_path / "b.tsv", "--samples", "A", "B",
            "--out-rel", tmp_path / "rel.tsv", "--out-counts", tmp_path / "cnt.tsv", "--out-qc", tmp_path / "qc.tsv")
    assert r.returncode == 0, r.stderr
    cnt = pd.read_csv(tmp_path / "cnt.tsv", sep="\t", index_col=0)
    assert cnt.loc["k__Bacteria|p__Firmicutes|c__C|o__O|f__F|g__Blautia|s__Blautia_obeum", "A"] == 800
    assert "UNCLASSIFIED" not in cnt.index
    assert pd.read_csv(tmp_path / "qc.tsv", sep="\t").set_index("sample_id").loc["A", "unclassified_pct"] == 20.0
    (tmp_path / "bad.tsv").write_text("#clade_name\ttax\trelative_abundance\nk__Bacteria\t2\t100\n")
    r = run("merge_metaphlan.py", "--profiles", tmp_path / "bad.tsv", "--samples", "X", "--out-rel", tmp_path / "r2",
            "--out-counts", tmp_path / "c2", "--out-qc", tmp_path / "q2")
    assert r.returncode != 0 and "not found" in r.stderr   # format drift is loud, not silent


def test_mag_pipeline_glue(tmp_path):
    # two samples with bins
    for s, bins in (("S1", [("metabat2.1", 95, 1.0), ("maxbin2.2", 60, 3.0), ("concoct.3", 40, 2.0)]),
                    ("S2", [("metabat2.1", 92, 3.0), ("metabat2.2", 70, 20.0)])):
        d = tmp_path / s / "bins"; d.mkdir(parents=True)
        rows = []
        for n, c, k in bins:
            (d / f"{n}.fa").write_text(f">{s}_c1\nACGT\n")
            rows.append(f"{n}\t{c}\t{k}\t0.9\t1000\t300\t2000000\t0.5\t1800\tNone")
        (tmp_path / s / "quality_report.tsv").write_text(
            "Name\tCompleteness\tContamination\tCoding_Density\tContig_N50\tAverage_Gene_Length\tGenome_Size\tGC_Content\tTotal_Coding_Sequences\tAdditional_Notes\n"
            + "\n".join(rows) + "\n")
    out = tmp_path / "coll"
    r = run("collect_bins.py", "--samples", "S1", "S2", "--bin-dirs", tmp_path / "S1/bins", tmp_path / "S2/bins",
            "--reports", tmp_path / "S1/quality_report.tsv", tmp_path / "S2/quality_report.tsv", "--outdir", out)
    assert r.returncode == 0, r.stderr
    q = pd.read_csv(out / "all_bins_quality.tsv", sep="\t").set_index("genome")
    assert q.loc["S1__metabat2.1.fa", "quality"] == "HQ" and q.loc["S1__maxbin2.2.fa", "quality"] == "MQ"
    assert q.loc["S1__concoct.3.fa", "quality"] == "low" and q.loc["S2__metabat2.2.fa", "quality"] == "low"
    assert sorted(p.name for p in (out / "passing").glob("*.fa")) == ["S1__maxbin2.2.fa", "S1__metabat2.1.fa", "S2__metabat2.1.fa"]
    assert len(pd.read_csv(out / "genome_info.csv")) == 3

    reps = tmp_path / "reps"; reps.mkdir()
    for n in ("S1__metabat2.1", "S1__maxbin2.2"):
        (reps / f"{n}.fa").write_text(">c\nACGT\n")
    gt = tmp_path / "gtdb"; gt.mkdir()
    (gt / "gtdbtk.bac120.summary.tsv").write_text(
        "user_genome\tclassification\tfastani_reference\tfastani_ani\tclassification_method\n"
        "S1__metabat2.1\td__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Blautia_A;s__Blautia_A obeum\tGCF_1\t98.2\tANI\n"
        "S1__maxbin2.2\td__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Blautia_A;s__\tN/A\tN/A\ttaxonomic novelty determined using RED\n")
    r = run("summarise_mags.py", "--rep-dir", reps, "--quality", out / "all_bins_quality.tsv", "--gtdbtk-dir", gt, "--out", tmp_path / "mags.tsv")
    assert r.returncode == 0, r.stderr
    m = pd.read_csv(tmp_path / "mags.tsv", sep="\t", na_values=["NA"], keep_default_na=False).set_index("mag")
    assert m.loc["S1__metabat2.1", "fastani_ani"] == pytest.approx(98.2) and pd.isna(m.loc["S1__maxbin2.2", "fastani_ani"])
    assert m.loc["S1__maxbin2.2", "gtdb_classification"].endswith("s__")

    rgi, abr = tmp_path / "rgi", tmp_path / "abr"; rgi.mkdir(); abr.mkdir()
    (rgi / "S1__metabat2.1.rgi.txt").write_text(
        "ORF_ID\tContig\tCut_Off\tBest_Hit_ARO\tBest_Identities\tDrug Class\n"
        "o1\tc1\tStrict\tacrF\t98.1\tfluoroquinolone antibiotic\no2\tc1\tLoose\tfoo\t40.0\tx\n")
    (abr / "S1__metabat2.1.ncbi.tsv").write_text(
        "#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n"
        "f\tc1\t1\t100\t+\tblaTEM-1\t1-100/100\t=====\t0/0\t100\t99\tncbi\tX\tbeta-lactamase\tBETA-LACTAM\n")
    r = run("summarise_amr.py", "--rgi-dir", rgi, "--abricate-dir", abr, "--out", tmp_path / "amr.tsv")
    assert r.returncode == 0, r.stderr
    a = pd.read_csv(tmp_path / "amr.tsv", sep="\t")
    assert sorted(a.gene) == ["acrF", "blaTEM-1"] and set(a.mag) == {"S1__metabat2.1"}


def test_merge_coverm_and_bracken(tmp_path):
    for s, vals in (("A", (60.0, 30.0, 10.0)), ("B", (50.0, 0.0, 50.0))):
        (tmp_path / f"{s}.cov").write_text(
            f"Genome\t{s}_R1 Relative Abundance (%)\t{s}_R1 Covered Fraction\nunmapped\t{vals[2]}\tNA\nMAG1\t{vals[0]}\t0.9\nMAG2\t{vals[1]}\t0.5\n")
    r = run("merge_coverm.py", "--inputs", tmp_path / "A.cov", tmp_path / "B.cov", "--samples", "A", "B",
            "--out-abundance", tmp_path / "ab.tsv", "--out-covered", tmp_path / "cv.tsv")
    assert r.returncode == 0, r.stderr
    ab = pd.read_csv(tmp_path / "ab.tsv", sep="\t", index_col=0)
    assert ab.loc["MAG1", "A"] == 60.0 and ab.loc["unmapped", "B"] == 50.0
    for s in "AB":
        (tmp_path / f"{s}.bracken").write_text("name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n"
                                               "Escherichia coli\t562\tS\t10\t5\t15\t0.5\nBlautia obeum\t40520\tS\t10\t5\t15\t0.5\n")
    r = run("merge_bracken.py", "--inputs", tmp_path / "A.bracken", tmp_path / "B.bracken", "--samples", "A", "B", "--out", tmp_path / "br.tsv")
    assert r.returncode == 0, r.stderr
    assert pd.read_csv(tmp_path / "br.tsv", sep="\t", index_col=0).loc["Escherichia coli", "A"] == 15


def test_build_sample_sheet_hmp2(tmp_path):
    meta = pd.DataFrame({
        "External ID": ["E1", "E1", "E2", "E3", "E4"], "Participant ID": ["P1", "P1", "P2", "P3", "P4"],
        "data_type": ["metagenomics", "16S", "metagenomics", "metagenomics", "viromics"],
        "week_num": ["0", "0", "4", "2", "1"], "diagnosis": ["CD", "CD", "UC", "weird", "UC"],
        "consent_age": ["30", "30", "40", "50", "20"], "Antibiotics": ["No"] * 5,
        "Immunosuppressants (e.g. oral corticosteroids)": ["Yes"] * 5, "hbi": ["3", "3", "NA", "1", "1"], "sccai": ["NA"] * 5})
    meta.to_csv(tmp_path / "meta.csv", index=False)
    ena = pd.DataFrame({
        "run_accession": ["R1", "R2", "R3", "R4"], "sample_alias": ["E1", "E1", "E2", "E9"],
        "library_strategy": ["WGS", "AMPLICON", "WGS", "WGS"], "instrument_model": ["HiSeq"] * 4,
        "first_created": ["2018-01-02"] * 4, "experiment_alias": ["x"] * 4,
        "fastq_ftp": ["h/a_1.fq.gz;h/a_2.fq.gz", "h/b_1.fq.gz;h/b_2.fq.gz", "h/c.fq.gz", "h/d_1.fq.gz;h/d_2.fq.gz"]})
    ena.to_csv(tmp_path / "ena.tsv", sep="\t", index=False)
    r = run("build_sample_sheet_hmp2.py", "--metadata", tmp_path / "meta.csv", "--ena", tmp_path / "ena.tsv",
            "--mapping", ROOT / "config" / "hmp2_columns.yaml", "--out", tmp_path / "samples.tsv")
    assert r.returncode == 0, r.stderr
    s = pd.read_csv(tmp_path / "samples.tsv", sep="\t")
    assert sorted(s.sample_id) == ["E1_16S", "E1_MGX"]                       # E2 single-end, E3 unmapped dx, E4 unmapped platform
    assert s.fastq_1.str.startswith("https://h/").all()
    sys.path.insert(0, str(ROOT / "workflow" / "lib"))
    import mbstats as mb
    assert len(mb.load_samples(tmp_path / "samples.tsv")) == 2              # output satisfies the workflow's own validator
    r = run("build_sample_sheet_hmp2.py", "--metadata", tmp_path / "meta.csv", "--ena", tmp_path / "ena.tsv",
            "--mapping", ROOT / "config" / "hmp2_columns.yaml", "--filter", "nonexistent=1", "--out", tmp_path / "x.tsv")
    assert r.returncode != 0 and "Available" in r.stderr


def test_method_agreement(tmp_path):
    import numpy as np
    rng = np.random.default_rng(0)
    sp = ["Escherichia coli", "Blautia obeum", "Bacteroides vulgatus", "Roseburia hominis", "Prevotella copri", "Akkermansia muciniphila"]
    base = rng.lognormal(3, 1, (len(sp), 12))
    a = pd.DataFrame(base * rng.lognormal(0, 0.05, base.shape), index=sp, columns=[f"S{i}" for i in range(12)])
    b = pd.DataFrame(base * rng.lognormal(0, 0.05, base.shape), index=[s.replace(" ", " ") for s in sp], columns=a.columns)
    a.index.name = b.index.name = "feature"
    a.to_csv(tmp_path / "mpa.tsv", sep="\t"); b.to_csv(tmp_path / "brk.tsv", sep="\t")
    r = run("method_agreement.py", "--metaphlan", tmp_path / "mpa.tsv", "--bracken", tmp_path / "brk.tsv", "--outdir", tmp_path / "o")
    assert r.returncode == 0, r.stderr
    s = pd.read_csv(tmp_path / "o" / "profiler_agreement_summary.tsv", sep="\t")
    assert s.n_shared_species.iloc[0] == 6 and s.median_species_rho.iloc[0] > 0.8

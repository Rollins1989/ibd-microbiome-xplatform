# Downstream analysis: metadata -> layers -> diversity / DA / concordance / longitudinal / networks / ML / MAGs -> figures.
import yaml

L = f"{OUT}/layers"
FIG = f"{OUT}/figures"
CFG = f"{OUT}/config/analysis.yaml"
R = f"{SCRIPTS}/R"
TAXA_LAYERS = [l for l in ACTIVE_LAYERS if config["layers"][l]["type"] == "counts"]

FINAL_TARGETS = [
    f"{FIG}/Fig1_study_design.pdf",
    f"{FIG}/Fig2_diversity.pdf",
    f"{FIG}/Fig3_platform_concordance.pdf",
    f"{FIG}/Fig4_differential_abundance.pdf",
    f"{FIG}/Fig5_ml_performance.pdf",
    f"{FIG}/Fig7_longitudinal.pdf",
    f"{FIG}/Fig8_cooccurrence_networks.pdf",
    f"{OUT}/report/summary.md",
]
if not config["analysis_only"]:
    FINAL_TARGETS += [f"{OUT}/shotgun/qc/multiqc_report.html"]
if USE_MAGS:
    FINAL_TARGETS += [f"{FIG}/Fig6_MAG_tree.pdf"]
if CLADES:
    FINAL_TARGETS += [f"{FIG}/Fig9_strain_persistence.pdf"]
if USE_KRAKEN:
    FINAL_TARGETS += [f"{OUT}/profiler_agreement/profiler_agreement_summary.tsv"]


rule analysis_config:
    """Snapshot of the configuration sections that influence the analysis (changing e.g. DADA2 settings does not rerun it)."""
    output:
        CFG,
    params:
        cfg={k: config[k] for k in ("design", "stats", "layers", "networks", "ml")}
        | {"amplicon": {"decontam": config["amplicon"]["decontam"]}, "mags": config["mags"]},
    run:
        with open(output[0], "w") as fh:
            yaml.safe_dump(params.cfg, fh, sort_keys=False)


rule prepare_metadata:
    input:
        samples=config["samples"],
        cfg=CFG,
    output:
        samples=f"{OUT}/metadata/samples_meta.tsv",
        specimens=f"{OUT}/metadata/specimens.tsv",
        events=f"{OUT}/metadata/flare_events.tsv",
    conda:
        f"{ENVS}/py_stats.yaml"
    log:
        f"{OUT}/logs/prepare_metadata.log",
    shell:
        """
        python {SCRIPTS}/prepare_metadata.py --samples {input.samples} --config {input.cfg} \
            --out-samples {output.samples} --out-specimens {output.specimens} --out-events {output.events} 2> {log}
        """


rule decontam:
    input:
        samples=config["samples"],
        asv=f"{OUT}/16s/dada2/asv_table.tsv",
        cfg=CFG,
    output:
        table=f"{OUT}/16s/asv_table.decontam.tsv",
        removed=f"{OUT}/16s/decontam_removed_asvs.tsv",
    conda:
        f"{ENVS}/r_stats.yaml"
    log:
        f"{OUT}/logs/decontam.log",
    shell:
        """
        Rscript {R}/decontam.R --samples {input.samples} --asv {input.asv} --config {input.cfg} \
            --out {output.table} --removed {output.removed} > {log} 2>&1
        """


_layer_outputs = {
    "meta": f"{L}/metadata.analysis.tsv",
    "specimens": f"{L}/specimens.analysis.tsv",
    "asv": f"{L}/asv_table.filtered.tsv",
    "summary": f"{L}/layer_summary.tsv",
    "dropped": f"{L}/dropped_samples.tsv",
}
for _l in ACTIVE_LAYERS:
    _layer_outputs[_l] = f"{L}/{_l}.tsv"
if USE_MAGS:
    _layer_outputs["mag"] = f"{L}/mag.tsv"


rule prepare_layers:
    input:
        meta=f"{OUT}/metadata/samples_meta.tsv",
        specimens=f"{OUT}/metadata/specimens.tsv",
        asv=f"{OUT}/16s/asv_table.decontam.tsv",
        tax=f"{OUT}/16s/dada2/taxonomy.tsv",
        mpa=f"{OUT}/shotgun/metaphlan/merged_estcounts.tsv",
        picrust=f"{OUT}/16s/picrust2/path_abun_unstrat.tsv",
        humann=([f"{OUT}/shotgun/humann/merged_pathabundance_cpm_unstratified.tsv"] if USE_HUMANN else []),
        mag=([f"{OUT}/mags/abundance/mag_relabund.tsv"] if USE_MAGS else []),
    output:
        **_layer_outputs,
    params:
        humann=(f"--humann {OUT}/shotgun/humann/merged_pathabundance_cpm_unstratified.tsv" if USE_HUMANN else ""),
        mag=(f"--mag-abundance {OUT}/mags/abundance/mag_relabund.tsv" if USE_MAGS else ""),
        min16=config["amplicon"]["min_reads_per_sample"],
        minmgx=config["shotgun"]["min_microbial_reads"],
        outdir=L,
    conda:
        f"{ENVS}/py_stats.yaml"
    log:
        f"{OUT}/logs/prepare_layers.log",
    shell:
        """
        python {SCRIPTS}/prepare_layers.py --samples-meta {input.meta} --specimens {input.specimens} --asv {input.asv} \
            --taxonomy {input.tax} --mpa-counts {input.mpa} --picrust {input.picrust} {params.humann} {params.mag} \
            --min-reads-16s {params.min16} --min-reads-mgx {params.minmgx} --out-dir {params.outdir} 2> {log}
        """


rule fig1_design:
    input:
        specimens=f"{L}/specimens.analysis.tsv",
        cfg=CFG,
    output:
        f"{FIG}/Fig1_study_design.pdf",
        f"{OUT}/design/design_summary.tsv",
    params:
        out=f"{OUT}/design",
        fig=FIG,
    conda:
        f"{ENVS}/r_stats.yaml"
    log:
        f"{OUT}/logs/fig1.log",
    shell:
        "Rscript {R}/fig1_design.R --config {input.cfg} --specimens {input.specimens} --figdir {params.fig} --outdir {params.out} > {log} 2>&1"


rule diversity:
    input:
        meta=f"{L}/metadata.analysis.tsv",
        asv=f"{L}/asv_table.filtered.tsv",
        mgx=f"{L}/mgx_species.tsv",
        tree=f"{OUT}/16s/tree/asv_tree.rooted.nwk",
        cfg=CFG,
    output:
        f"{FIG}/Fig2_diversity.pdf",
        f"{OUT}/diversity/alpha_diversity.tsv",
        f"{OUT}/diversity/alpha_models.tsv",
        f"{OUT}/diversity/permanova_pairwise.tsv",
    params:
        layers=L,
        out=f"{OUT}/diversity",
        fig=FIG,
    threads: config["stats"]["threads"]
    conda:
        f"{ENVS}/r_stats.yaml"
    log:
        f"{OUT}/logs/diversity.log",
    shell:
        """
        Rscript {R}/diversity.R --config {input.cfg} --meta {input.meta} --layer-dir {params.layers} --tree {input.tree} \
            --outdir {params.out} --figdir {params.fig} > {log} 2>&1
        """


def _da_inputs(wc):
    return {"layer": f"{L}/{wc.layer}.tsv", "meta": f"{L}/metadata.analysis.tsv", "cfg": CFG}


_DA_ENV = {"clr_lmm": "r_stats", "maaslin2": "da_maaslin2", "ancombc2": "da_ancombc2", "aldex2": "da_aldex2"}
for _m, _env in _DA_ENV.items():

    rule:
        name:
            f"da_{_m}"
        input:
            unpack(_da_inputs),
        output:
            f"{OUT}/da/{{layer}}.{_m}.tsv",
        params:
            type=lambda wc: config["layers"][wc.layer]["type"],
            script=f"{R}/da_{_m}.R",
        threads: config["stats"]["threads"]
        conda:
            f"{ENVS}/{_env}.yaml"
        log:
            f"{OUT}/logs/da/{{layer}}.{_m}.log",
        shell:
            """
            Rscript {params.script} --config {input.cfg} --meta {input.meta} --layer {input.layer} \
                --layer-name {wildcards.layer} --type {params.type} --out {output} --threads {threads} > {log} 2>&1
            """


rule da_consensus:
    input:
        lambda wc: [f"{OUT}/da/{wc.layer}.{m}.tsv" for m in config["layers"][wc.layer]["da_methods"]],
    output:
        f"{OUT}/da/{{layer}}.consensus.tsv",
    params:
        q=config["stats"]["q_threshold"],
        k=config["stats"]["min_methods"],
    conda:
        f"{ENVS}/py_stats.yaml"
    shell:
        "python {SCRIPTS}/da_consensus.py --inputs {input} --q {params.q} --min-methods {params.k} --out {output}"


rule fig4_da_heatmap:
    input:
        consensus=expand(f"{OUT}/da/{{l}}.consensus.tsv", l=ACTIVE_LAYERS),
        cfg=CFG,
    output:
        f"{FIG}/Fig4_differential_abundance.pdf",
    params:
        inputs=lambda wc, input: ",".join(input.consensus),
        fig=FIG,
    conda:
        f"{ENVS}/r_stats.yaml"
    log:
        f"{OUT}/logs/fig4.log",
    shell:
        "Rscript {R}/fig_da_heatmap.R --config {input.cfg} --inputs {params.inputs} --figdir {params.fig} > {log} 2>&1"


rule concordance:
    input:
        specimens=f"{L}/specimens.analysis.tsv",
        meta=f"{L}/metadata.analysis.tsv",
        layers=expand(f"{L}/{{l}}.tsv", l=ACTIVE_LAYERS),
        cfg=CFG,
    output:
        f"{FIG}/Fig3_platform_concordance.pdf",
        f"{OUT}/concordance/genus_spearman.tsv",
    params:
        layers=L,
        out=f"{OUT}/concordance",
        fig=FIG,
    conda:
        f"{ENVS}/r_stats.yaml"
    log:
        f"{OUT}/logs/concordance.log",
    shell:
        """
        Rscript {R}/concordance.R --config {input.cfg} --specimens {input.specimens} --meta {input.meta} \
            --layer-dir {params.layers} --outdir {params.out} --figdir {params.fig} > {log} 2>&1
        """


rule longitudinal:
    input:
        meta=f"{L}/metadata.analysis.tsv",
        specimens=f"{L}/specimens.analysis.tsv",
        events=f"{OUT}/metadata/flare_events.tsv",
        alpha=f"{OUT}/diversity/alpha_diversity.tsv",
        layers=[f"{L}/16s_genus.tsv", f"{L}/mgx_species.tsv"],
        cfg=CFG,
    output:
        f"{FIG}/Fig7_longitudinal.pdf",
        f"{OUT}/longitudinal/activity_associations.tsv",
    params:
        layers=L,
        out=f"{OUT}/longitudinal",
        fig=FIG,
    conda:
        f"{ENVS}/r_stats.yaml"
    log:
        f"{OUT}/logs/longitudinal.log",
    shell:
        """
        Rscript {R}/longitudinal.R --config {input.cfg} --meta {input.meta} --specimens {input.specimens} \
            --events {input.events} --alpha {input.alpha} --layer-dir {params.layers} \
            --outdir {params.out} --figdir {params.fig} > {log} 2>&1
        """


rule networks:
    input:
        meta=f"{L}/metadata.analysis.tsv",
        layer=f"{L}/{config['networks']['layer']}.tsv",
        cfg=CFG,
    output:
        f"{FIG}/Fig8_cooccurrence_networks.pdf",
        f"{OUT}/networks/network_metrics.tsv",
    params:
        layers=L,
        out=f"{OUT}/networks",
        fig=FIG,
    conda:
        f"{ENVS}/py_stats.yaml"
    log:
        f"{OUT}/logs/networks.log",
    shell:
        "python {SCRIPTS}/network.py --config {input.cfg} --meta {input.meta} --layer-dir {params.layers} --outdir {params.out} --figdir {params.fig} 2> {log}"


rule ml:
    input:
        specimens=f"{L}/specimens.analysis.tsv",
        layers=expand(f"{L}/{{l}}.tsv", l=["16s_genus", "mgx_species"]),
        cfg=CFG,
    output:
        f"{FIG}/Fig5_ml_performance.pdf",
        f"{OUT}/ml/ml_summary.tsv",
    params:
        layers=L,
        out=f"{OUT}/ml",
        fig=FIG,
    threads: config["ml"]["n_jobs"]
    conda:
        f"{ENVS}/py_stats.yaml"
    log:
        f"{OUT}/logs/ml.log",
    shell:
        "python {SCRIPTS}/ml.py --config {input.cfg} --specimens {input.specimens} --layer-dir {params.layers} --outdir {params.out} --figdir {params.fig} 2> {log}"


rule ml_permutation_control:
    """Negative control: AUROC on subject-shuffled labels must hover around 0.5 (run on demand)."""
    input:
        specimens=f"{L}/specimens.analysis.tsv",
        layers=expand(f"{L}/{{l}}.tsv", l=["16s_genus", "mgx_species"]),
        cfg=CFG,
    output:
        f"{OUT}/ml_permuted/ml_summary.tsv",
    params:
        layers=L,
        out=f"{OUT}/ml_permuted",
    conda:
        f"{ENVS}/py_stats.yaml"
    shell:
        "python {SCRIPTS}/ml.py --permute-labels --config {input.cfg} --specimens {input.specimens} --layer-dir {params.layers} --outdir {params.out} --figdir {params.out}/fig"


rule mag_analysis:
    input:
        meta=f"{L}/metadata.analysis.tsv",
        mag=f"{L}/mag.tsv",
        mpa=f"{L}/mgx_species.tsv",
        summary=f"{OUT}/mags/summary/mag_summary.tsv",
        amr=f"{OUT}/mags/summary/amr_hits.tsv",
        tree=f"{OUT}/mags/tree/mag_tree.nwk",
        cfg=CFG,
    output:
        f"{FIG}/Fig6_MAG_tree.pdf",
        f"{OUT}/mags/analysis/mag_catalogue_annotated.tsv",
        f"{OUT}/mags/analysis/mag_summary_stats.tsv",
    params:
        layers=L,
        out=f"{OUT}/mags/analysis",
        fig=FIG,
    threads: config["stats"]["threads"]
    conda:
        f"{ENVS}/r_stats.yaml"
    log:
        f"{OUT}/logs/mag_analysis.log",
    shell:
        """
        Rscript {R}/mag_analysis.R --config {input.cfg} --meta {input.meta} --layer-dir {params.layers} \
            --mag-summary {input.summary} --amr {input.amr} --tree {input.tree} --outdir {params.out} \
            --figdir {params.fig} --threads {threads} > {log} 2>&1
        """


rule strain_analysis:
    input:
        meta=f"{L}/metadata.analysis.tsv",
        trees=expand(f"{OUT}/shotgun/strainphlan/{{c}}/tree.nwk", c=CLADES),
    output:
        f"{FIG}/Fig9_strain_persistence.pdf",
        f"{OUT}/strain/strain_summary.tsv",
    params:
        out=f"{OUT}/strain",
        fig=FIG,
    conda:
        f"{ENVS}/py_stats.yaml"
    shell:
        "python {SCRIPTS}/strain_analysis.py --meta {input.meta} --trees {input.trees} --outdir {params.out} --figdir {params.fig}"


rule profiler_agreement:
    input:
        mpa=f"{L}/mgx_species.tsv",
        brk=f"{OUT}/shotgun/bracken/merged_species_counts.tsv",
    output:
        f"{OUT}/profiler_agreement/profiler_agreement_summary.tsv",
    params:
        out=f"{OUT}/profiler_agreement",
    conda:
        f"{ENVS}/py_stats.yaml"
    shell:
        "python {SCRIPTS}/method_agreement.py --metaphlan {input.mpa} --bracken {input.brk} --outdir {params.out}"


rule summary_report:
    input:
        [t for t in FINAL_TARGETS if t.endswith(".pdf")],
        consensus=expand(f"{OUT}/da/{{l}}.consensus.tsv", l=ACTIVE_LAYERS),
    output:
        f"{OUT}/report/summary.md",
    params:
        out=OUT,
    conda:
        f"{ENVS}/py_stats.yaml"
    shell:
        "python {SCRIPTS}/make_report.py --outdir {params.out} --out {output}"

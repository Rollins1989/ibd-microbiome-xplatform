# Shotgun branch: MetaPhlAn 4 -> HUMAnN -> StrainPhlAn; optional Kraken2/Bracken.
MPA_DB, MPA_IDX = config["databases"]["metaphlan_dir"], config["databases"]["metaphlan_index"]


rule metaphlan:
    input:
        r1=f"{OUT}/shotgun/clean/{{sample}}_R1.fastq.gz",
        r2=f"{OUT}/shotgun/clean/{{sample}}_R2.fastq.gz",
    output:
        stats=f"{OUT}/shotgun/metaphlan/{{sample}}.profile_stats.tsv",
        profile=f"{OUT}/shotgun/metaphlan/{{sample}}.profile.tsv",
        bt2=f"{OUT}/shotgun/metaphlan/{{sample}}.bowtie2.bz2",
        sam=(f"{OUT}/shotgun/metaphlan/{{sample}}.sam.bz2" if CLADES else temp(f"{OUT}/shotgun/metaphlan/{{sample}}.sam.bz2")),
    wildcard_constraints:
        sample=MGX_RE,
    params:
        db=MPA_DB,
        idx=MPA_IDX,
        sam_flag=lambda wc, output: f"-s {output.sam}" if CLADES else "",
        touch_sam=":" if CLADES else "touch",
    threads: 16
    resources:
        mem_mb=24000,
    conda:
        f"{ENVS}/biobakery.yaml"
    log:
        f"{OUT}/logs/metaphlan/{{sample}}.log",
    shell:
        """
        # run 1: mapping + read statistics (counts are needed for count-based DA methods)
        metaphlan {input.r1},{input.r2} --input_type fastq --bowtie2db {params.db} --index {params.idx} \
            --nproc {threads} -t rel_ab_w_read_stats --bowtie2out {output.bt2} {params.sam_flag} \
            -o {output.stats} > {log} 2>&1
        # run 2 (re-uses the bowtie2 output, takes seconds): standard profile required by HUMAnN
        metaphlan {output.bt2} --input_type bowtie2out --bowtie2db {params.db} --index {params.idx} \
            --nproc {threads} -t rel_ab -o {output.profile} >> {log} 2>&1
        {params.touch_sam} {output.sam}
        """


rule merge_metaphlan:
    input:
        expand(f"{OUT}/shotgun/metaphlan/{{s}}.profile_stats.tsv", s=MGX),
    output:
        rel=f"{OUT}/shotgun/metaphlan/merged_relabund.tsv",
        counts=f"{OUT}/shotgun/metaphlan/merged_estcounts.tsv",
        qc=f"{OUT}/shotgun/metaphlan/unclassified_qc.tsv",
    params:
        samples=" ".join(MGX),
    conda:
        f"{ENVS}/py_stats.yaml"
    log:
        f"{OUT}/logs/merge_metaphlan.log",
    shell:
        """
        python {SCRIPTS}/merge_metaphlan.py --profiles {input} --samples {params.samples} \
            --out-rel {output.rel} --out-counts {output.counts} --out-qc {output.qc} 2> {log}
        """


rule humann:
    input:
        r1=f"{OUT}/shotgun/clean/{{sample}}_R1.fastq.gz",
        r2=f"{OUT}/shotgun/clean/{{sample}}_R2.fastq.gz",
        profile=f"{OUT}/shotgun/metaphlan/{{sample}}.profile.tsv",
    output:
        genes=f"{OUT}/shotgun/humann/samples/{{sample}}_genefamilies.tsv",
        paths=f"{OUT}/shotgun/humann/samples/{{sample}}_pathabundance.tsv",
        cover=f"{OUT}/shotgun/humann/samples/{{sample}}_pathcoverage.tsv",
    wildcard_constraints:
        sample=MGX_RE,
    params:
        outdir=f"{OUT}/shotgun/humann/samples",
        nuc=config["databases"]["humann_nucleotide"],
        prot=config["databases"]["humann_protein"],
    threads: 16
    resources:
        mem_mb=48000,
    conda:
        f"{ENVS}/biobakery.yaml"
    log:
        f"{OUT}/logs/humann/{{sample}}.log",
    shell:
        """
        tmp=$(mktemp -d); trap 'rm -rf $tmp' EXIT
        cat {input.r1} {input.r2} > $tmp/{wildcards.sample}.fastq.gz     # concatenated gzip members are valid gzip
        humann --input $tmp/{wildcards.sample}.fastq.gz --output {params.outdir} \
               --output-basename {wildcards.sample} --taxonomic-profile {input.profile} \
               --nucleotide-database {params.nuc} --protein-database {params.prot} \
               --threads {threads} --remove-temp-output > {log} 2>&1
        """


rule humann_merge:
    input:
        paths=expand(f"{OUT}/shotgun/humann/samples/{{s}}_pathabundance.tsv", s=MGX),
        genes=expand(f"{OUT}/shotgun/humann/samples/{{s}}_genefamilies.tsv", s=MGX),
    output:
        unstrat=f"{OUT}/shotgun/humann/merged_pathabundance_cpm_unstratified.tsv",
        strat=f"{OUT}/shotgun/humann/merged_pathabundance_cpm_stratified.tsv",
        genes=f"{OUT}/shotgun/humann/merged_genefamilies_cpm.tsv",
    params:
        indir=f"{OUT}/shotgun/humann/samples",
        outdir=f"{OUT}/shotgun/humann",
    conda:
        f"{ENVS}/biobakery.yaml"
    log:
        f"{OUT}/logs/humann_merge.log",
    shell:
        """
        humann_join_tables --input {params.indir} --output {params.outdir}/merged_pathabundance.tsv --file_name pathabundance > {log} 2>&1
        humann_renorm_table --input {params.outdir}/merged_pathabundance.tsv --output {params.outdir}/merged_pathabundance_cpm.tsv \
            --units cpm --update-snames >> {log} 2>&1
        humann_split_stratified_table --input {params.outdir}/merged_pathabundance_cpm.tsv --output {params.outdir} >> {log} 2>&1
        humann_join_tables --input {params.indir} --output {params.outdir}/merged_genefamilies.tsv --file_name genefamilies >> {log} 2>&1
        humann_renorm_table --input {params.outdir}/merged_genefamilies.tsv --output {output.genes} --units cpm --update-snames >> {log} 2>&1
        """


# ------------------------------------------------------------------ StrainPhlAn (optional)
rule sample2markers:
    input:
        f"{OUT}/shotgun/metaphlan/{{sample}}.sam.bz2",
    output:
        f"{OUT}/shotgun/strainphlan/markers/{{sample}}.pkl",
    wildcard_constraints:
        sample=MGX_RE,
    params:
        db=(f"-d {config['shotgun']['strainphlan']['database']}" if config["shotgun"]["strainphlan"]["database"] else ""),
    threads: 4
    conda:
        f"{ENVS}/biobakery.yaml"
    log:
        f"{OUT}/logs/strainphlan/{{sample}}.sample2markers.log",
    shell:
        """
        tmp=$(mktemp -d); trap 'rm -rf $tmp' EXIT
        sample2markers.py -i {input} -o $tmp -n {threads} {params.db} > {log} 2>&1
        mv $tmp/*.pkl {output}      # tip label = sample id (file name)
        """


rule strainphlan:
    input:
        expand(f"{OUT}/shotgun/strainphlan/markers/{{s}}.pkl", s=MGX),
    output:
        tree=f"{OUT}/shotgun/strainphlan/{{clade}}/tree.nwk",
    params:
        d=f"{OUT}/shotgun/strainphlan/{{clade}}",
        n=config["shotgun"]["strainphlan"]["marker_in_n_samples"],
        m=config["shotgun"]["strainphlan"]["sample_with_n_markers"],
        db=(f"-d {config['shotgun']['strainphlan']['database']}" if config["shotgun"]["strainphlan"]["database"] else ""),
    threads: 16
    conda:
        f"{ENVS}/biobakery.yaml"
    log:
        f"{OUT}/logs/strainphlan/{{clade}}.log",
    shell:
        """
        strainphlan -s {input} -o {params.d} -c {wildcards.clade} -n {threads} \
            --marker_in_n_samples {params.n} --sample_with_n_markers {params.m} \
            --phylophlan_mode accurate {params.db} > {log} 2>&1
        cp "$(ls {params.d}/RAxML_bestTree.*.tre | head -n 1)" {output.tree}
        """


# ------------------------------------------------------------------ Kraken2 / Bracken (optional)
rule kraken2:
    input:
        r1=f"{OUT}/shotgun/clean/{{sample}}_R1.fastq.gz",
        r2=f"{OUT}/shotgun/clean/{{sample}}_R2.fastq.gz",
    output:
        report=f"{OUT}/shotgun/kraken2/{{sample}}.kreport",
    wildcard_constraints:
        sample=MGX_RE,
    params:
        db=config["databases"]["kraken2"],
    threads: 16
    resources:
        mem_mb=96000,
    conda:
        f"{ENVS}/kraken.yaml"
    log:
        f"{OUT}/logs/kraken2/{{sample}}.log",
    shell:
        """
        kraken2 --db {params.db} --paired --gzip-compressed --threads {threads} --report {output.report} \
                --output /dev/null {input.r1} {input.r2} > {log} 2>&1
        """


rule bracken:
    input:
        f"{OUT}/shotgun/kraken2/{{sample}}.kreport",
    output:
        f"{OUT}/shotgun/bracken/{{sample}}.bracken",
    wildcard_constraints:
        sample=MGX_RE,
    params:
        db=config["databases"]["kraken2"],
        r=config["shotgun"]["kraken2"]["read_length"],
        lvl=config["shotgun"]["kraken2"]["bracken_level"],
        thr=config["shotgun"]["kraken2"]["bracken_threshold"],
        rep=lambda wc: f"{OUT}/shotgun/bracken/{wc.sample}.bracken_report",
    conda:
        f"{ENVS}/kraken.yaml"
    log:
        f"{OUT}/logs/bracken/{{sample}}.log",
    shell:
        "bracken -d {params.db} -i {input} -o {output} -w {params.rep} -r {params.r} -l {params.lvl} -t {params.thr} > {log} 2>&1"


rule merge_bracken:
    input:
        expand(f"{OUT}/shotgun/bracken/{{s}}.bracken", s=MGX),
    output:
        f"{OUT}/shotgun/bracken/merged_species_counts.tsv",
    params:
        samples=" ".join(MGX),
    conda:
        f"{ENVS}/py_stats.yaml"
    shell:
        "python {SCRIPTS}/merge_bracken.py --inputs {input} --samples {params.samples} --out {output}"

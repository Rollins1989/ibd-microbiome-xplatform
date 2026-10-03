# 16S branch: primer removal -> DADA2 (per batch) -> taxonomy -> phylogeny -> PICRUSt2.


def s16_input(sample):
    d = "16s/trimmed" if PRIMERS else "raw"
    return [f"{OUT}/{d}/{sample}_R1.fastq.gz", f"{OUT}/{d}/{sample}_R2.fastq.gz"]


def batch_samples(batch):
    return [s for s in S16_ALL if SIDX.loc[s, "batch"] == batch]


rule cutadapt_16s:
    input:
        unpack(raw_reads),
    output:
        r1=f"{OUT}/16s/trimmed/{{sample}}_R1.fastq.gz",
        r2=f"{OUT}/16s/trimmed/{{sample}}_R2.fastq.gz",
        report=f"{OUT}/16s/trimmed/{{sample}}.cutadapt.json",
    wildcard_constraints:
        sample=S16_RE,
    params:
        fwd=config["amplicon"]["primer_fwd"],
        rev=config["amplicon"]["primer_rev"],
    threads: 4
    conda:
        f"{ENVS}/qc.yaml"
    log:
        f"{OUT}/logs/cutadapt/{{sample}}.log",
    shell:
        """
        cutadapt -g {params.fwd} -G {params.rev} --discard-untrimmed --minimum-length 50 -j {threads} \
                 --json {output.report} -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log} 2>&1
        """


rule dada2_batch:
    input:
        lambda wc: [f for s in batch_samples(wc.batch) for f in s16_input(s)],
    output:
        seqtab=f"{OUT}/16s/dada2/batches/{{batch}}.seqtab.rds",
        track=f"{OUT}/16s/dada2/batches/{{batch}}.track.tsv",
        errplot=f"{OUT}/16s/dada2/batches/{{batch}}.error_models.pdf",
    wildcard_constraints:
        batch=BATCH_RE,
    params:
        ids=lambda wc: ",".join(batch_samples(wc.batch)),
        r1=lambda wc: ",".join(s16_input(s)[0] for s in batch_samples(wc.batch)),
        r2=lambda wc: ",".join(s16_input(s)[1] for s in batch_samples(wc.batch)),
        filtdir=lambda wc: f"{OUT}/16s/dada2/filtered/{wc.batch}",
        tl=",".join(map(str, config["amplicon"]["dada2"]["trunc_len"])),
        ee=",".join(map(str, config["amplicon"]["dada2"]["max_ee"])),
        tq=config["amplicon"]["dada2"]["trunc_q"],
        ov=config["amplicon"]["dada2"]["min_overlap"],
        mm=config["amplicon"]["dada2"]["max_mismatch"],
        pool=str(config["amplicon"]["dada2"]["pool"]).lower(),
        seed=config["stats"]["seed"],
    threads: 16
    resources:
        mem_mb=32000,
    conda:
        f"{ENVS}/dada2.yaml"
    log:
        f"{OUT}/logs/dada2/{{batch}}.log",
    shell:
        """
        Rscript {SCRIPTS}/R/dada2_batch.R --samples {params.ids} --r1 {params.r1} --r2 {params.r2} \
            --filtdir {params.filtdir} --trunc-len {params.tl} --max-ee {params.ee} --trunc-q {params.tq} \
            --min-overlap {params.ov} --max-mismatch {params.mm} --pool {params.pool} --threads {threads} \
            --seed {params.seed} --out-seqtab {output.seqtab} --out-track {output.track} \
            --out-errplot {output.errplot} > {log} 2>&1
        rm -rf {params.filtdir}
        """


rule dada2_merge:
    input:
        seqtabs=expand(f"{OUT}/16s/dada2/batches/{{b}}.seqtab.rds", b=BATCHES_16S),
        tracks=expand(f"{OUT}/16s/dada2/batches/{{b}}.track.tsv", b=BATCHES_16S),
    output:
        asv=f"{OUT}/16s/dada2/asv_table.tsv",
        fasta=f"{OUT}/16s/dada2/asv_seqs.fasta",
        tax=f"{OUT}/16s/dada2/taxonomy.tsv",
        track=f"{OUT}/16s/dada2/read_tracking.tsv",
    params:
        seqtabs=lambda wc, input: ",".join(input.seqtabs),
        tracks=lambda wc, input: ",".join(input.tracks),
        train=config["databases"]["dada2_train_set"],
        species=config["databases"]["dada2_species_set"],
        boot=config["amplicon"]["dada2"]["min_boot"],
        minlen=config["amplicon"]["dada2"]["min_asv_length"],
        outdir=f"{OUT}/16s/dada2",
        seed=config["stats"]["seed"],
    threads: 16
    resources:
        mem_mb=32000,
    conda:
        f"{ENVS}/dada2.yaml"
    log:
        f"{OUT}/logs/dada2/merge.log",
    shell:
        """
        Rscript {SCRIPTS}/R/dada2_merge.R --seqtabs {params.seqtabs} --tracks {params.tracks} \
            --train-set {params.train} --species-set "{params.species}" --min-boot {params.boot} \
            --min-length {params.minlen} --threads {threads} --seed {params.seed} --outdir {params.outdir} > {log} 2>&1
        """


rule filter_fasta:
    """Restrict the ASV FASTA to ASVs that survived decontamination."""
    input:
        fasta=f"{OUT}/16s/dada2/asv_seqs.fasta",
        table=f"{OUT}/16s/asv_table.decontam.tsv",
    output:
        f"{OUT}/16s/asv_seqs.decontam.fasta",
    run:
        keep = set(l.split("\t", 1)[0] for l in open(input.table).read().splitlines()[1:])
        write, n = False, 0
        with open(input.fasta) as fi, open(output[0], "w") as fo:
            for line in fi:
                if line.startswith(">"):
                    write = line[1:].strip().split()[0] in keep
                    n += write
                if write:
                    fo.write(line)
        assert n == len(keep), f"FASTA lacks {len(keep) - n} ASVs present in the table"


rule mafft:
    input:
        f"{OUT}/16s/asv_seqs.decontam.fasta",
    output:
        f"{OUT}/16s/tree/asv_aligned.fasta",
    threads: 16
    conda:
        f"{ENVS}/phylogeny.yaml"
    log:
        f"{OUT}/logs/mafft.log",
    shell:
        "mafft --auto --thread {threads} {input} > {output} 2> {log}"


rule fasttree:
    input:
        f"{OUT}/16s/tree/asv_aligned.fasta",
    output:
        f"{OUT}/16s/tree/asv_tree.unrooted.nwk",
    conda:
        f"{ENVS}/phylogeny.yaml"
    log:
        f"{OUT}/logs/fasttree.log",
    shell:
        "FastTree -nt -gtr -gamma {input} > {output} 2> {log}"


rule root_tree:
    input:
        f"{OUT}/16s/tree/asv_tree.unrooted.nwk",
    output:
        f"{OUT}/16s/tree/asv_tree.rooted.nwk",
    conda:
        f"{ENVS}/r_stats.yaml"
    log:
        f"{OUT}/logs/root_tree.log",
    shell:
        "Rscript {SCRIPTS}/R/root_tree.R --tree {input} --out {output} > {log} 2>&1"


rule picrust2:
    input:
        fasta=f"{OUT}/16s/asv_seqs.decontam.fasta",
        table=f"{OUT}/16s/asv_table.decontam.tsv",
    output:
        f"{OUT}/16s/picrust2/path_abun_unstrat.tsv",
    params:
        run=f"{OUT}/16s/picrust2/run",
    threads: 16
    conda:
        f"{ENVS}/picrust2.yaml"
    log:
        f"{OUT}/logs/picrust2.log",
    shell:
        """
        rm -rf {params.run}   # PICRUSt2 refuses to write into an existing directory
        picrust2_pipeline.py -s {input.fasta} -i {input.table} -o {params.run} -p {threads} > {log} 2>&1
        gunzip -c {params.run}/pathways_out/path_abun_unstrat.tsv.gz > {output}
        """

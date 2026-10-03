# Raw read retrieval and shotgun QC / host depletion.

rule get_reads:
    """Fetch (URL) or link (local path) the raw FASTQ named in the sample sheet."""
    output:
        f"{OUT}/raw/{{sample}}_R{{mate}}.fastq.gz",
    params:
        src=lambda wc: str(SIDX.loc[wc.sample, f"fastq_{wc.mate}"]),
    log:
        f"{OUT}/logs/get_reads/{{sample}}_R{{mate}}.log",
    threads: 1
    retries: 2
    shell:
        r"""
        if [[ "{params.src}" =~ ^(https?|ftp):// ]]; then
            curl -fL --retry 5 --retry-delay 10 -o {output} "{params.src}" 2> {log}
            gzip -t {output}
        else
            ln -sf "$(realpath "{params.src}")" {output}
        fi
        """


rule fastp_mgx:
    input:
        unpack(raw_reads),
    output:
        r1=temp(f"{OUT}/shotgun/qc/{{sample}}_R1.fastq.gz"),
        r2=temp(f"{OUT}/shotgun/qc/{{sample}}_R2.fastq.gz"),
        json=f"{OUT}/shotgun/qc/{{sample}}.fastp.json",
        html=f"{OUT}/shotgun/qc/{{sample}}.fastp.html",
    wildcard_constraints:
        sample=MGX_RE,
    params:
        minlen=config["shotgun"]["fastp"]["min_length"],
        q=config["shotgun"]["fastp"]["cut_mean_quality"],
    threads: 8
    resources:
        mem_mb=8000,
    conda:
        f"{ENVS}/qc.yaml"
    log:
        f"{OUT}/logs/fastp/{{sample}}.log",
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
              --detect_adapter_for_pe --cut_tail --cut_window_size 4 --cut_mean_quality {params.q} \
              --length_required {params.minlen} --n_base_limit 5 --low_complexity_filter \
              --thread {threads} --json {output.json} --html {output.html} 2> {log}
        """


rule host_removal:
    """Map to GRCh38 and keep read pairs where both mates are unmapped."""
    input:
        r1=f"{OUT}/shotgun/qc/{{sample}}_R1.fastq.gz",
        r2=f"{OUT}/shotgun/qc/{{sample}}_R2.fastq.gz",
    output:
        r1=f"{OUT}/shotgun/clean/{{sample}}_R1.fastq.gz",
        r2=f"{OUT}/shotgun/clean/{{sample}}_R2.fastq.gz",
    wildcard_constraints:
        sample=MGX_RE,
    params:
        idx=config["databases"]["host_bowtie2_index"],
    threads: 16
    resources:
        mem_mb=16000,
    conda:
        f"{ENVS}/mapping.yaml"
    log:
        bt2=f"{OUT}/logs/host_removal/{{sample}}.bowtie2.log",
    shell:
        """
        bowtie2 -p {threads} --very-sensitive-local -x {params.idx} -1 {input.r1} -2 {input.r2} 2> {log.bt2} \
          | samtools view -@ 2 -b -f 12 -F 256 - \
          | samtools fastq -@ 2 -n -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -
        """


rule multiqc_shotgun:
    input:
        expand(f"{OUT}/shotgun/qc/{{s}}.fastp.json", s=MGX),
        expand(f"{OUT}/logs/host_removal/{{s}}.bowtie2.log", s=MGX),
    output:
        f"{OUT}/shotgun/qc/multiqc_report.html",
    conda:
        f"{ENVS}/qc.yaml"
    log:
        f"{OUT}/logs/multiqc_shotgun.log",
    shell:
        "multiqc --force -o $(dirname {output}) -n multiqc_report.html {OUT}/shotgun/qc {OUT}/logs/host_removal > {log} 2>&1"

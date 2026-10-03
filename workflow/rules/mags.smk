# Metagenome-assembled genomes: assembly -> 3 binners -> DAS Tool -> CheckM2 -> dRep -> GTDB-Tk -> annotation/AMR -> abundance.
MC = config["mags"]
M = f"{OUT}/mags"


rule megahit:
    input:
        r1=f"{OUT}/shotgun/clean/{{sample}}_R1.fastq.gz",
        r2=f"{OUT}/shotgun/clean/{{sample}}_R2.fastq.gz",
    output:
        f"{M}/assembly/{{sample}}/contigs.fa",
    wildcard_constraints:
        sample=ASM_RE,
    params:
        tmp=lambda wc: f"{M}/assembly/{wc.sample}.tmp",
        minlen=MC["min_contig_len"],
    threads: 32
    resources:
        mem_mb=128000,
    conda:
        f"{ENVS}/assembly.yaml"
    log:
        f"{OUT}/logs/megahit/{{sample}}.log",
    shell:
        r"""
        rm -rf {params.tmp}      # MEGAHIT refuses an existing output directory
        megahit -1 {input.r1} -2 {input.r2} -t {threads} --min-contig-len {params.minlen} -o {params.tmp} > {log} 2>&1
        # unique, space-free contig names across samples
        sed -E 's/^>([^ ]+).*/>{wildcards.sample}_\1/' {params.tmp}/final.contigs.fa > {output}
        rm -rf {params.tmp}
        """


rule map_for_binning:
    input:
        contigs=f"{M}/assembly/{{sample}}/contigs.fa",
        r1=f"{OUT}/shotgun/clean/{{sample}}_R1.fastq.gz",
        r2=f"{OUT}/shotgun/clean/{{sample}}_R2.fastq.gz",
    output:
        bam=f"{M}/binning/{{sample}}/reads.sorted.bam",
        bai=f"{M}/binning/{{sample}}/reads.sorted.bam.bai",
        depth=f"{M}/binning/{{sample}}/depth.txt",
    wildcard_constraints:
        sample=ASM_RE,
    params:
        idx=lambda wc: f"{M}/binning/{wc.sample}/contigs_index",
    threads: 16
    resources:
        mem_mb=32000,
    conda:
        f"{ENVS}/binning.yaml"
    log:
        f"{OUT}/logs/binning/{{sample}}.map.log",
    shell:
        """
        bowtie2-build --threads {threads} {input.contigs} {params.idx} > {log} 2>&1
        bowtie2 -p {threads} --no-unal -x {params.idx} -1 {input.r1} -2 {input.r2} 2>> {log} \
            | samtools sort -@ 4 -o {output.bam} -
        samtools index {output.bam}
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {output.bam} >> {log} 2>&1
        rm -f {params.idx}*.bt2*
        """


rule metabat2:
    input:
        contigs=f"{M}/assembly/{{sample}}/contigs.fa",
        depth=f"{M}/binning/{{sample}}/depth.txt",
    output:
        directory(f"{M}/binning/{{sample}}/metabat2"),
    wildcard_constraints:
        sample=ASM_RE,
    params:
        minlen=max(1500, MC["min_contig_len"]),
    threads: 8
    conda:
        f"{ENVS}/binning.yaml"
    log:
        f"{OUT}/logs/binning/{{sample}}.metabat2.log",
    shell:
        """
        mkdir -p {output}
        metabat2 -i {input.contigs} -a {input.depth} -o {output}/bin -t {threads} -m {params.minlen} --seed 1 > {log} 2>&1 \
            || echo "WARNING: MetaBAT2 returned non-zero (typically: no bins formed)" >> {log}
        """


rule maxbin2:
    input:
        contigs=f"{M}/assembly/{{sample}}/contigs.fa",
        depth=f"{M}/binning/{{sample}}/depth.txt",
    output:
        directory(f"{M}/binning/{{sample}}/maxbin2"),
    wildcard_constraints:
        sample=ASM_RE,
    threads: 8
    conda:
        f"{ENVS}/binning.yaml"
    log:
        f"{OUT}/logs/binning/{{sample}}.maxbin2.log",
    shell:
        """
        mkdir -p {output}
        tail -n +2 {input.depth} | cut -f1,3 > {output}/abundance.txt
        run_MaxBin.pl -contig {input.contigs} -abund {output}/abundance.txt -out {output}/bin -thread {threads} > {log} 2>&1 \
            || echo "WARNING: MaxBin2 returned non-zero (typically: no bins formed)" >> {log}
        """


rule concoct:
    input:
        contigs=f"{M}/assembly/{{sample}}/contigs.fa",
        bam=f"{M}/binning/{{sample}}/reads.sorted.bam",
    output:
        directory(f"{M}/binning/{{sample}}/concoct"),
    wildcard_constraints:
        sample=ASM_RE,
    threads: 8
    conda:
        f"{ENVS}/binning.yaml"
    log:
        f"{OUT}/logs/binning/{{sample}}.concoct.log",
    shell:
        """
        mkdir -p {output}/work {output}/bins
        cut_up_fasta.py {input.contigs} -c 10000 -o 0 --merge_last -b {output}/work/contigs_10K.bed > {output}/work/contigs_10K.fa 2> {log}
        concoct_coverage_table.py {output}/work/contigs_10K.bed {input.bam} > {output}/work/coverage.tsv 2>> {log}
        concoct --composition_file {output}/work/contigs_10K.fa --coverage_file {output}/work/coverage.tsv \
                -b {output}/work/ -t {threads} >> {log} 2>&1 \
            || echo "WARNING: CONCOCT returned non-zero" >> {log}
        if [ -s {output}/work/clustering_gt1000.csv ]; then
            merge_cutup_clustering.py {output}/work/clustering_gt1000.csv > {output}/work/clustering_merged.csv 2>> {log}
            extract_fasta_bins.py {input.contigs} {output}/work/clustering_merged.csv --output_path {output}/bins >> {log} 2>&1
        fi
        """


rule dastool:
    input:
        contigs=f"{M}/assembly/{{sample}}/contigs.fa",
        mb=f"{M}/binning/{{sample}}/metabat2",
        mx=f"{M}/binning/{{sample}}/maxbin2",
        cc=f"{M}/binning/{{sample}}/concoct",
    output:
        directory(f"{M}/dastool/{{sample}}/bins"),
    wildcard_constraints:
        sample=ASM_RE,
    params:
        d=lambda wc: f"{M}/dastool/{wc.sample}",
        score=0.5,
    threads: 8
    conda:
        f"{ENVS}/binning.yaml"
    log:
        f"{OUT}/logs/dastool/{{sample}}.log",
    shell:
        r"""
        d={params.d}; mkdir -p $d; : > {log}
        Fasta_to_Contigs2Bin.sh -i {input.mb}      -e fa    > $d/metabat2.tsv 2>> {log} || true
        Fasta_to_Contigs2Bin.sh -i {input.mx}      -e fasta > $d/maxbin2.tsv  2>> {log} || true
        Fasta_to_Contigs2Bin.sh -i {input.cc}/bins -e fa    > $d/concoct.tsv  2>> {log} || true
        lst=""; names=""
        for b in metabat2 maxbin2 concoct; do
            if [ -s $d/$b.tsv ]; then lst="$lst,$d/$b.tsv"; names="$names,$b"; fi
        done
        if [ -n "$lst" ]; then
            DAS_Tool -i ${{lst#,}} -l ${{names#,}} -c {input.contigs} -o $d/DASTool --write_bins \
                     --score_threshold {params.score} --threads {threads} >> {log} 2>&1 \
                || echo "WARNING: DAS_Tool returned non-zero (typically: no bin above score threshold)" >> {log}
        else
            echo "WARNING: no binner produced bins for {wildcards.sample}" >> {log}
        fi
        if [ -d $d/DASTool_DASTool_bins ]; then mv $d/DASTool_DASTool_bins {output}; else mkdir -p {output}; fi
        """


rule checkm2:
    input:
        f"{M}/dastool/{{sample}}/bins",
    output:
        f"{M}/checkm2/{{sample}}/quality_report.tsv",
    wildcard_constraints:
        sample=ASM_RE,
    params:
        db=config["databases"]["checkm2"],
        out=lambda wc: f"{M}/checkm2/{wc.sample}",
    threads: 16
    resources:
        mem_mb=32000,
    conda:
        f"{ENVS}/checkm2.yaml"
    log:
        f"{OUT}/logs/checkm2/{{sample}}.log",
    shell:
        r"""
        if ls {input}/*.fa > /dev/null 2>&1; then
            checkm2 predict -i {input} -x fa -o {params.out} -t {threads} --database_path {params.db} --force > {log} 2>&1
        else
            echo "no bins for {wildcards.sample}" > {log}
            printf 'Name\tCompleteness\tContamination\tCoding_Density\tContig_N50\tAverage_Gene_Length\tGenome_Size\tGC_Content\tTotal_Coding_Sequences\tAdditional_Notes\n' > {output}
        fi
        """


rule collect_bins:
    input:
        dirs=expand(f"{M}/dastool/{{s}}/bins", s=ASM),
        reports=expand(f"{M}/checkm2/{{s}}/quality_report.tsv", s=ASM),
    output:
        passing=directory(f"{M}/bins_passing/passing"),
        info=f"{M}/bins_passing/genome_info.csv",
        quality=f"{M}/bins_passing/all_bins_quality.tsv",
    params:
        samples=" ".join(ASM),
        outdir=f"{M}/bins_passing",
        mc=MC["min_completeness"],
        mx=MC["max_contamination"],
        hc=MC["hq_completeness"],
        hx=MC["hq_contamination"],
    conda:
        f"{ENVS}/py_stats.yaml"
    shell:
        """
        python {SCRIPTS}/collect_bins.py --samples {params.samples} --bin-dirs {input.dirs} --reports {input.reports} \
            --min-completeness {params.mc} --max-contamination {params.mx} --hq-completeness {params.hc} \
            --hq-contamination {params.hx} --outdir {params.outdir}
        """


rule drep:
    input:
        passing=f"{M}/bins_passing/passing",
        info=f"{M}/bins_passing/genome_info.csv",
    output:
        directory(f"{M}/drep/dereplicated_genomes"),
    params:
        out=f"{M}/drep",
        ani=MC["derep_ani"],
        mc=MC["min_completeness"],
        mx=MC["max_contamination"],
    threads: 16
    resources:
        mem_mb=32000,
    conda:
        f"{ENVS}/mag_tools.yaml"
    log:
        f"{OUT}/logs/drep.log",
    shell:
        """
        n=$(ls {input.passing}/*.fa 2>/dev/null | wc -l)
        if [ "$n" -lt 2 ]; then echo "dRep needs >= 2 passing MAGs, found $n" >&2; exit 1; fi
        dRep dereplicate {params.out} -g {input.passing}/*.fa -p {threads} -pa 0.9 -sa {params.ani} \
             --genomeInfo {input.info} -comp {params.mc} -con {params.mx} > {log} 2>&1
        """


rule gtdbtk:
    input:
        f"{M}/drep/dereplicated_genomes",
    output:
        touch(f"{M}/gtdbtk/.classify_done"),
    params:
        out=f"{M}/gtdbtk",
        db=config["databases"]["gtdbtk"],
    threads: 32
    resources:
        mem_mb=128000,
    conda:
        f"{ENVS}/gtdbtk.yaml"
    log:
        f"{OUT}/logs/gtdbtk_classify.log",
    shell:
        """
        export GTDBTK_DATA_PATH={params.db}
        gtdbtk classify_wf --genome_dir {input} --out_dir {params.out} -x fa --cpus {threads} --skip_ani_screen > {log} 2>&1
        ls {params.out}/gtdbtk.*.summary.tsv > /dev/null
        """


rule mag_tree:
    """De novo tree of the representative MAGs (bacterial 120 markers) from the GTDB-Tk alignment, midpoint rooted."""
    input:
        reps=f"{M}/drep/dereplicated_genomes",
        done=f"{M}/gtdbtk/.classify_done",
    output:
        f"{M}/tree/mag_tree.nwk",
    params:
        d=f"{M}/tree/work",
        db=config["databases"]["gtdbtk"],
    threads: 16
    conda:
        f"{ENVS}/gtdbtk.yaml"
    log:
        f"{OUT}/logs/mag_tree.log",
    shell:
        r"""
        export GTDBTK_DATA_PATH={params.db}
        rm -rf {params.d}; mkdir -p {params.d}
        gtdbtk identify --genome_dir {input.reps} --out_dir {params.d}/identify -x fa --cpus {threads} > {log} 2>&1
        gtdbtk align --identify_dir {params.d}/identify --out_dir {params.d}/align --skip_gtdb_refs --cpus {threads} >> {log} 2>&1
        msa=$(find {params.d}/align -name 'gtdbtk.bac120.user_msa.fasta*' | head -n 1)
        [ -n "$msa" ] || {{ echo "no bacterial MSA produced (archaea-only MAG set?)" >&2; exit 1; }}
        if [[ "$msa" == *.gz ]]; then zcat "$msa" > {params.d}/msa.faa; else cp "$msa" {params.d}/msa.faa; fi
        FastTree -lg -gamma {params.d}/msa.faa > {params.d}/unrooted.nwk 2>> {log}
        Rscript {SCRIPTS}/R/root_tree.R --tree {params.d}/unrooted.nwk --out {output} >> {log} 2>&1
        """


rule bakta:
    input:
        f"{M}/drep/dereplicated_genomes",
    output:
        touch(f"{M}/annotation/.bakta_done"),
    params:
        out=f"{M}/annotation",
        db=config["databases"]["bakta"],
    threads: 16
    conda:
        f"{ENVS}/annotation.yaml"
    log:
        f"{OUT}/logs/bakta.log",
    shell:
        """
        : > {log}
        for fa in {input}/*.fa; do
            n=$(basename $fa .fa)
            bakta --db {params.db} --output {params.out}/$n --prefix $n --threads {threads} --force $fa >> {log} 2>&1
        done
        """


rule amr_screen:
    input:
        f"{M}/drep/dereplicated_genomes",
    output:
        touch(f"{M}/amr/.screen_done"),
    params:
        out=f"{M}/amr",
    threads: 16
    conda:
        f"{ENVS}/amr.yaml"
    log:
        f"{OUT}/logs/amr_screen.log",
    shell:
        """
        mkdir -p {params.out}/rgi {params.out}/abricate; : > {log}
        for fa in {input}/*.fa; do
            n=$(basename $fa .fa)
            rgi main --input_sequence $fa --output_file {params.out}/rgi/$n.rgi -t contig -a DIAMOND -n {threads} --clean --low_quality >> {log} 2>&1
            for db in ncbi vfdb; do
                abricate --db $db --minid 80 --mincov 80 --threads {threads} $fa > {params.out}/abricate/$n.$db.tsv 2>> {log}
            done
        done
        """


rule mag_summary:
    input:
        reps=f"{M}/drep/dereplicated_genomes",
        quality=f"{M}/bins_passing/all_bins_quality.tsv",
        done=f"{M}/gtdbtk/.classify_done",
        amr=f"{M}/amr/.screen_done",
        bakta=f"{M}/annotation/.bakta_done",
    output:
        mags=f"{M}/summary/mag_summary.tsv",
        amr=f"{M}/summary/amr_hits.tsv",
    params:
        gtdb=f"{M}/gtdbtk",
        rgi=f"{M}/amr/rgi",
        abr=f"{M}/amr/abricate",
    conda:
        f"{ENVS}/py_stats.yaml"
    shell:
        """
        python {SCRIPTS}/summarise_mags.py --rep-dir {input.reps} --quality {input.quality} --gtdbtk-dir {params.gtdb} --out {output.mags}
        python {SCRIPTS}/summarise_amr.py --rgi-dir {params.rgi} --abricate-dir {params.abr} --out {output.amr}
        """


rule coverm:
    input:
        reps=f"{M}/drep/dereplicated_genomes",
        r1=f"{OUT}/shotgun/clean/{{sample}}_R1.fastq.gz",
        r2=f"{OUT}/shotgun/clean/{{sample}}_R2.fastq.gz",
    output:
        f"{M}/abundance/per_sample/{{sample}}.tsv",
    wildcard_constraints:
        sample=MGX_RE,
    params:
        ident=MC["coverm"]["min_identity"],
        aligned=MC["coverm"]["min_aligned"],
        cov=MC["coverm"]["min_covered_percent"],
    threads: 8
    resources:
        mem_mb=24000,
    conda:
        f"{ENVS}/mag_tools.yaml"
    log:
        f"{OUT}/logs/coverm/{{sample}}.log",
    shell:
        """
        coverm genome --coupled {input.r1} {input.r2} --genome-fasta-directory {input.reps} -x fa \
            --methods relative_abundance covered_fraction --min-read-percent-identity {params.ident} \
            --min-read-aligned-percent {params.aligned} --min-covered-fraction {params.cov} \
            -t {threads} -o {output} > {log} 2>&1
        """


rule merge_coverm:
    input:
        expand(f"{M}/abundance/per_sample/{{s}}.tsv", s=MGX),
    output:
        ab=f"{M}/abundance/mag_relabund.tsv",
        cov=f"{M}/abundance/mag_covered_fraction.tsv",
    params:
        samples=" ".join(MGX),
    conda:
        f"{ENVS}/py_stats.yaml"
    shell:
        "python {SCRIPTS}/merge_coverm.py --inputs {input} --samples {params.samples} --out-abundance {output.ab} --out-covered {output.cov}"

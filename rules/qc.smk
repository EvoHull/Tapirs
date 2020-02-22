# ==================================================
# QUALITY CONTROL SNAKEFILE
# ==================================================
# --------------------------------------------------
# fastp, control for sequence quality and pair reads
# --------------------------------------------------

configfile: "config.yaml"

rule fastp_trim_and_merge:
    message:
        "Beginning fastp quality control of raw data"
    conda:
        "../envs/tapirs.yaml"
    input:
        read1 = "data/01_demultiplexed/{library}/{sample}.R1.fastq.gz",
        read2 = "data/01_demultiplexed/{library}/{sample}.R2.fastq.gz"
    output:
        out1 = "results/02_trimmed/{library}/{sample}.R1.fastq.gz",
        out2 = "results/02_trimmed/{library}/{sample}.R2.fastq.gz",
        out_unpaired1 = "results/02_trimmed/{library}/{sample}.unpaired.R1.fastq.gz",
        out_unpaired2 = "results/02_trimmed/{library}/{sample}.unpaired.R2.fastq.gz",
        out_failed = "results/02_trimmed/{library}/{sample}.failed.fastq.gz",
        merged = "results/02_trimmed/{library}/{sample}.merged.fastq.gz",
        json = "reports/fastp/{library}/{sample}.json",
        html = "reports/fastp/{library}/{sample}.html"
    shell:
        "fastp \
        -i {input.read1} \
        -I {input.read2} \
        -o {output.out1} \
        -O {output.out2} \
        --unpaired1 {output.out_unpaired1} \
        --unpaired2 {output.out_unpaired2} \
        --failed_out {output.out_failed} \
        -j {output.json} \
        -h {output.html} \
        --qualified_quality_phred 30 \
        --length_required 90 \
        --cut_tail \
        --trim_front1 20 \
        --trim_front2 20 \
        --max_len1 106 \
        --max_len2 106 \
        --merge \
        --merged_out {output.merged} \
        --overlap_len_require 90 \
        --correction \
        "

# unpaired files need bringing back in

    #     cat 3_qc/${library}/${sample}.R1.fastq >> 4_merged_fastq/${library}/${sample}_merged.fastq\n",
    # "        cat 3_qc/${library}/${sample}_unpaired.R1.fastq >> 4_merged_fastq/${library}/${sample}_merged.fastq

# -----------------------------------------------------
# convert files from fastq to fasta
# -----------------------------------------------------

rule fastq_to_fasta:
    conda:
        "../envs/tapirs.yaml"
    input:
        "results/02_trimmed/{library}/{sample}.merged.fastq.gz"
    output:
        "results/02_trimmed/{library}/{sample}.merged.fasta",
    shell:
        "vsearch \
        --fastq_filter {input} \
        --fastaout {output} \
        "


# -----------------------------------------------------
# vsearch, fastq report
# -----------------------------------------------------

rule vsearch_fastq_report:
    conda:
        "../envs/tapirs.yaml"
    input:
        "results/02_trimmed/{library}/{sample}.merged.fastq.gz"
    output:
        fqreport = "reports/vsearch/{library}/{sample}_fq_eestats",
        fqreadstats = "reports/vsearch/{library}/{sample}_fq_readstats"
    shell:
        "vsearch \
        --fastq_eestats {input} \
        --output {output.fqreport} ; \
        vsearch \
        --fastq_stats {input} \
        --log {output.fqreadstats} \
        "


# -----------------------------------------------------
# dereplication
# -----------------------------------------------------

rule vsearch_dereplication:
    conda:
        "../envs/tapirs.yaml"
    input:
        "results/02_trimmed/{library}/{sample}.merged.fasta"
    output:
        "results/02_trimmed/{library}/{sample}.merged.derep.fasta"
    shell:
        "vsearch \
        --derep_fulllength {input} \
        --sizeout \
        --minuniquesize 3 \
        --output {output} \
        "


# -----------------------------------------------------
# denoise
# -----------------------------------------------------

rule vsearch_denoising:
    conda:
        "../envs/tapirs.yaml"
    input:
        "results/02_trimmed/{library}/{sample}.merged.derep.fasta"
    output:
        fasta = "results/03_denoised/{library}/{sample}.fasta",
        biom = "reports/vsearch/{library}/{sample}.denoise.biom"
    #params:
    #    log="reports/denoise/{library}/vsearch.log"
    shell:
        "vsearch \
        --cluster_unoise {input} \
        --centroids {output.fasta} \
        --biomout {output.biom} \
        --minsize 3 \
        --unoise_alpha 0.5 \
        "
        #" --notrunclabels
        # --log {params.log} \


# -----------------------------------------------------
# chimera removal
# -----------------------------------------------------

rule vsearch_dechimerisation: # output needs fixing
    conda:
        "../envs/tapirs.yaml"
    input:
        "results/03_denoised/{library}/{sample}.fasta"
    output: # fix
        text = "results/03_denoised/{library}/{sample}_chimera.txt",
        fasta = "results/03_denoised/{library}/nc_{sample}.fasta"
    shell:
        "vsearch \
        --uchime_ref {input} \
        --db data/databases/12S_full/12s_full.fasta \
        --mindiffs 1 \
        --mindiv 0.8 \
        --uchimeout {output.text} \
        --nonchimeras {output.fasta} \
        "


# ------------------------------------------------------
# re-replication
# -------------------------------------------------------

rule vsearch_rereplication:
    conda:
        "../envs/tapirs.yaml"
    input:
        "results/03_denoised/{library}/{sample}.fasta" # check
    output:
        "results/rereplicated/{library}/{sample}.fasta"
    threads:
        6
    shell:
        "vsearch \
        --rereplicate {input} \
        --output {output} \
        "

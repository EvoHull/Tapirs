# ==================================================
# QUALITY CONTROL SNAKEFILE
# ==================================================

configfile: "config.yaml"

# --------------------------------------------------
# fastp, trim on length and sequence quality
# --------------------------------------------------

rule fastp_trim_reads:
    conda:
        "../envs/environment.yaml"
    input:
        # read1 = "data/01_demultiplexed/test1/1EB.R1.fastq.gz",
        # read2 = "data/01_demultiplexed/test1/1EB.R2.fastq.gz",
        read1 = expand(
            "data/01_demultiplexed/{library}/{sample}.R1.fastq.gz", sample=SAMPLES, library=LIBRARIES),
        read2 = expand(
            "data/01_demultiplexed/{library}/{sample}.R2.fastq.gz", sample=SAMPLES, library=LIBRARIES),
    output:
        R1trimmed = "results/02_trimmed/{library}/{sample}.R1.trimmed.fastq",
        R2trimmed = "results/02_trimmed/{library}/{sample}.R2.trimmed.fastq",
        unpairedR1 = "results/02_trimmed/{library}/{sample}.R1.unpaired.fastq",
        unpairedR2 = "results/02_trimmed/{library}/{sample}.R2.unpaired.fastq",
        failed = "results/02_trimmed/{library}/{sample}.failedfilter.fastq",
        json = "reports/fastp/{library}/{sample}.fastp.json",  # report
        html = "reports/fastp/{library}/{sample}.fastp.html",  # report
    shell:
        "fastp \
        --in1 {input.read1} \
        --in2 {input.read2} \
        --out1 {output.R1trimmed} \
        --out2 {output.R2trimmed} \
        --unpaired1 {output.unpairedR1} \
        --unpaired2 {output.unpairedR2} \
        --failed_out {output.failed} \
        -j {output.json} \
        -h {output.html} \
        --qualified_quality_phred {config[FASTP_qual_phred]} \
        --length_required {config[FASTP_len_required]} \
        --cut_tail \
        --trim_front1 {config[FASTP_trim_front1]} \
        --trim_front2 {config[FASTP_trim_front2]} \
        --max_len1 {config[FASTP_max_len1]} \
        --max_len2 {config[FASTP_max_len2]} \
        --correction \
        -w 6"

# --------------------------------------------------
# fastp, merge reads R1 and R2
# merged.fastq also includes unmerged and unpaired
# --------------------------------------------------
rule fastp_merge_reads:
    conda:
        "../envs/environment.yaml"
    input:
        trimmedread1 = expand(
            "results/02_trimmed/{library}/{sample}.R1.trimmed.fastq", sample=SAMPLES, library=LIBRARIES),
        trimmedread2 = expand(
            "results/02_trimmed/{library}/{sample}.R2.trimmed.fastq", sample=SAMPLES, library=LIBRARIES),
        unpairedR1 = expand(
            "results/02_trimmed/{library}/{sample}.R1.unpaired.fastq", sample=SAMPLES, library=LIBRARIES),
        unpairedR2 = expand(
            "results/02_trimmed/{library}/{sample}.R2.unpaired.fastq", sample=SAMPLES, library=LIBRARIES),
    output:
        merged = "results/03_merged/{library}/{sample}.concat.fastq",
        # unmerged1 = "results/03_merged/{library}/{sample}.unmerged1.fastq",
        # unmerged2 = "results/03_merged/{library}/{sample}.unmerged2.fastq",
    params:
        overlap = config["FASTP_overlap_len"]
    shell:
        "fastp \
        --in1 {input.trimmedread1} \
        --in2 {input.trimmedread2} \
        --unpaired1 {input.unpairedR1} \
        --unpaired2 {input.unpairedR2} \
        --merge \
        --include_unmerged \
        --merged_out {output.merged} \
        --overlap_len_require {params.overlap} \
        "

# -----------------------------------------------------
# keep forward unpaired and convert fastq to fasta
#  not needed with fastp --keep-unmerged command?
# -----------------------------------------------------

# rule seqkit_merge_in_unpaired_to_fasta:
#     conda:
#         "../envs/environment.yaml"
#     input:
#         merged = "results/03_merged/{library}/{sample}.merged.fastq",
#         unpairedR1 = "results/03_merged/{library}/{sample}.unpaired.R1.fastq",
#         R1 = "results/02_trimmed/{library}/{sample}.R1.fastq"
#     output:
#         fa = "results/03_merged/{library}/{sample}.concat.fasta",
#         fq = "results/03_merged/{library}/{sample}.concat.fastq"
#     shell:
#         """
#         seqkit fq2fa \
#         {input.merged} {input.unpairedR1} {input.R1} \
#         -o {output.fa} {output.fq}
#         """

# -----------------------------------------------------
# find and remove empty merged files
#  -will work if files in /library directory ie -depth 2
# -----------------------------------------------------

rule remove_empty_files:
    input:
        "results/03_merged/"
    output:
        "reports/empty_files_deleted.txt"
    shell:
        "find {input} -depth 2 -type f -empty -print > {output} -delete"

# -----------------------------------------------------
# seqkit fastq to fasta
# -----------------------------------------------------

rule seqkit_convert_to_fasta:
    conda:
        "../envs/environment.yaml"
    input:
        expand("results/03_merged/{library}/{sample}.concat.fastq",
               sample=SAMPLES, library=LIBRARIES)
    output:
        "results/03_merged/{library}/{sample}.concat.fasta",
    shell:
        "seqkit fq2fa {input} -o {output}"

# -----------------------------------------------------
# vsearch dereplication
# -----------------------------------------------------

rule vsearch_dereplication:
    conda:
        "../envs/environment.yaml"
    input:
        expand("results/03_merged/{library}/{sample}.concat.fasta", sample=SAMPLES, library=LIBRARIES)
    output:
        "results/04_dereplicated/{library}/{sample}.derep.fasta"
    params:
        minuniqsize = config["VSEARCH_minuniqsize"]
    shell:
        "vsearch --derep_fulllength {input} --sizeout --minuniquesize {params.minuniqsize} --output {output}"

# -----------------------------------------------------
# vsearch denoise
# -----------------------------------------------------

rule vsearch_denoising:
    conda:
        "../envs/environment.yaml"
    input:
        expand(
            "results/04_dereplicated/{library}/{sample}.derep.fasta", sample=SAMPLES, library=LIBRARIES),
    output:
        centroids = "results/05_denoised/{library}/{sample}.denoise.fasta",
        cluster_results = "reports/vsearch/{library}/{sample}.denoise-report.txt"
    shell:
        "vsearch \
        --cluster_unoise {input} \
        --sizein \
        --sizeout \
        --minsize {config[VSEARCH_minsize]} \
        --unoise_alpha {config[VSEARCH_unoise_alpha]} \
        --id {config[VSEARCH_id]} \
        --centroids {output.centroids} \
        --uc {output.cluster_results}"

# -----------------------------------------------------
# vsearch chimera removal
# -----------------------------------------------------

rule vsearch_dechimerisation:
    conda:
        "../envs/environment.yaml"
    input:
        seqs = expand(
            "results/05_denoised/{library}/{sample}.denoise.fasta", sample=SAMPLES, library=LIBRARIES),
        blast_db = config["dechim_blast_db"]
    output:
        chimeras = "results/06_dechimera/{library}/{sample}.chimera.fasta.gz",
        borderline = "results/06_dechimera/{library}/{sample}.borderline-chimera.fasta.gz",
        nonchimeras = "results/06_dechimera/{library}/{sample}_nonchimera.fasta"
    params:
        mindiffs = config["VSEARCH_mindiffs"],
        mindiv = config["VSEARCH_mindiv"]
    shell:
        "vsearch \
        --uchime_ref {input.seqs} \
        --db {input.blast_db} \
        --mindiffs {params.mindiffs} \
        --mindiv {params.mindiv} \
        --chimeras {output.chimeras} \
        --borderline {output.borderline} \
        --nonchimeras {output.nonchimeras} \
        "

# ------------------------------------------------------
# vsearch re-replication
# -------------------------------------------------------

rule vsearch_rereplication:
    conda:
        "../envs/environment.yaml"
    input:
        expand("results/06_dechimera/{library}/{sample}_nonchimera.fasta",
               sample=SAMPLES, library=LIBRARIES)
    output:
        "results/07_rereplicated/{library}/{sample}.rerep.fasta"
    threads:
        6
    shell:
        "vsearch --rereplicate {input} --output {output}"

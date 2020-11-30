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
        read1 = "data/01_demultiplexed/{SAMPLES}.R1.fastq.gz",
        read2 = "data/01_demultiplexed/{SAMPLES}.R2.fastq.gz"
    output:
        R1trimmed = "results/02_trimmed/{SAMPLES}.R1.trimmed.fastq",
        R2trimmed = "results/02_trimmed/{SAMPLES}.R2.trimmed.fastq",
        R1unpaired = "results/02_trimmed/{SAMPLES}.R1.unpaired.fastq",
        R2unpaired = "results/02_trimmed/{SAMPLES}.R2.unpaired.fastq",
        failed = "results/02_trimmed/{SAMPLES}.failed.fastq",
        json = "reports/fastp/{SAMPLES}.fastp.json",
        html = "reports/fastp/{SAMPLES}.fastp.html"
    threads:
        10
    shell:
         "fastp \
        --in1 {input.read1} \
        --in2 {input.read2} \
        --out1 {output.R1trimmed} \
        --out2 {output.R2trimmed} \
        --unpaired1 {output.R1unpaired} \
        --unpaired2 {output.R2unpaired} \
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
        -w {config[threads]}"

# ------------------------------------------------------------------------------
# FASTP MERGE PAIRED END READS

rule fastp_merge_reads:
    conda:
        "../envs/environment.yaml"
    input:
        R1trimmed = "results/02_trimmed/{SAMPLES}.R1.trimmed.fastq",
        R2trimmed = "results/02_trimmed/{SAMPLES}.R2.trimmed.fastq",
    output:
        merged = "results/03_merged/{SAMPLES}.merged.fastq",
        R1unmerged = "results/03_merged/{SAMPLES}.R1.unmerged.fastq",
        R2unmerged = "results/03_merged/{SAMPLES}.R2.unmerged.fastq",
        json = "reports/fastp_merged/{SAMPLES}.merged.fastp.json",
        html = "reports/fastp_merged/{SAMPLES}.merged.fastp.html"
    shell:
        "fastp \
        --disable_quality_filtering \
        --disable_adapter_trimming \
        --in1 {input.R1trimmed} \
        --in2 {input.R2trimmed} \
        --out1 {output.R1unmerged} \
        --out2 {output.R2unmerged} \
        --merge \
        --merged_out {output.merged} \
        --overlap_len_require {config[FASTP_overlap_len]} \
        -w {config[threads]} \
        -j {output.json} \
        -h {output.html} "

# ------------------------------------------------------------------------------
# MERGE FORWARD READS

rule merge_forward_reads:
    # conda:
    #     "../envs/environment.yaml"
    input:
        merged = "results/03_merged/{SAMPLES}.merged.fastq",
        R1unpaired = "results/02_trimmed/{SAMPLES}.R1.unpaired.fastq",
        R1unmerged = "results/03_merged/{SAMPLES}.R1.unmerged.fastq"
    output:
        fq = "results/04_forward_merged/{SAMPLES}.forward.merged.fastq"
    run:
        filenames = [input.merged, input.R1unpaired, input.R1unmerged]
        with open(str(output.fq), 'w') as outfile:
            for fname in filenames:
                with open(str(fname)) as infile:
                    outfile.write(infile.read())

# ------------------------------------------------------------------------------
# FASTQ TO FASTA

rule seqkit_fq2fa:
    conda:
        "../envs/environment.yaml"
    input:
        fq = "results/04_forward_merged/{SAMPLES}.forward.merged.fastq"
    output:
        fa = "results/05_forward_merged/{SAMPLES}.fasta"
    shell:
        "seqkit fq2fa \
        {input.fq} \
        -o {output.fa}"

# ------------------------------------------------------------------------------

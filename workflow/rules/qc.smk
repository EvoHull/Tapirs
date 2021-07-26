# ==================================================
# QUALITY CONTROL SNAKEFILE
# ==================================================

# --------------------------------------------------
# fastp, trim on length and sequence quality
# --------------------------------------------------

rule fastp_trim_reads:
    conda:
        config['conda']
    input:
        read1 = config['input_data'] + "/{LIBRARIES}/{SAMPLES}.R1.fastq.gz",
        read2 = config['input_data'] + "/{LIBRARIES}/{SAMPLES}.R2.fastq.gz"
    params:
        overlap = int(config['FASTP_max_len1'] * ((config['FASTP_min_overlap']) / 100)),
    output:
        R1trimmed = "results/02_trimmed/{LIBRARIES}/{SAMPLES}.R1.trimmed.fastq",
        R2trimmed = "results/02_trimmed/{LIBRARIES}/{SAMPLES}.R2.trimmed.fastq",
        R1unpaired = "results/02_trimmed/{LIBRARIES}/{SAMPLES}.R1.unpaired.fastq",
        R2unpaired = "results/02_trimmed/{LIBRARIES}/{SAMPLES}.R2.unpaired.fastq",
        failed = "results/02_trimmed/{LIBRARIES}/{SAMPLES}.trimmed.failed.fastq",
        json = "results/02_trimmed/fastp_trimmed_reports/{LIBRARIES}/{SAMPLES}.fastp.json",
        html = "results/02_trimmed/fastp_trimmed_reports/{LIBRARIES}/{SAMPLES}.fastp.html"
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
        --unqualified_percent_limit {config[Fastp_unqualified_percent_limit]} \
        --average_qual {config[FASTP_qual_phred]} \
        --cut_front \
        --cut_tail \
        --cut_window_size {config[FASTP_window_size]} \
        --length_required {config[FASTP_len_required]} \
        --trim_front1 {config[FASTP_trim_front1]} \
        --trim_front2 {config[FASTP_trim_front2]} \
        --overlap_len_require {params.overlap} \
        --overlap_diff_percent_limit {config[FASTP_max_diff]} \
        --max_len1 {config[FASTP_max_len1]} \
        --max_len2 {config[FASTP_max_len2]}"

# ------------------------------------------------------------------------------
# FASTP MERGE PAIRED END READS

rule fastp_merge_reads:
    conda:
        config['conda']
    input:
        R1trimmed = "results/02_trimmed/{LIBRARIES}/{SAMPLES}.R1.trimmed.fastq",
        R2trimmed = "results/02_trimmed/{LIBRARIES}/{SAMPLES}.R2.trimmed.fastq",
    params:
        overlap = int(config['FASTP_max_len1'] * ((config['FASTP_min_overlap']) / 100)),
    output:
        merged = "results/03_merged/{LIBRARIES}/{SAMPLES}.merged.fastq",
        R1unmerged = "results/03_merged/{LIBRARIES}/{SAMPLES}.R1.unmerged.fastq",
        R2unmerged = "results/03_merged/{LIBRARIES}/{SAMPLES}.R2.unmerged.fastq",
        json = "results/03_merged/fastp_merged_reports/{LIBRARIES}/{SAMPLES}.merged.fastp.json",
        html = "results/03_merged/fastp_merged_reports/{LIBRARIES}/{SAMPLES}.merged.fastp.html"
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
        --overlap_len_require {params.overlap} \
        --overlap_diff_percent_limit {config[FASTP_max_diff]} \
        --length_limit {config[FASTP_length_limit]} \
        -j {output.json} \
        -h {output.html} \
        --correction"

# ------------------------------------------------------------------------------
# MERGE FORWARD READS

rule merge_forward_reads:
    input:
        merged = "results/03_merged/{LIBRARIES}/{SAMPLES}.merged.fastq",
        R1unpaired = "results/02_trimmed/{LIBRARIES}/{SAMPLES}.R1.unpaired.fastq",
        R1unmerged = "results/03_merged/{LIBRARIES}/{SAMPLES}.R1.unmerged.fastq"
    output:
        fq = "results/04_forward_merged/{LIBRARIES}/{SAMPLES}.forward.merged.fastq"
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
        config['conda']
    input:
        fq = "results/04_forward_merged/{LIBRARIES}/{SAMPLES}.forward.merged.fastq"
    output:
        fa = "results/05_forward_merged/{LIBRARIES}/{SAMPLES}.fasta"
    shell:
        "seqkit fq2fa \
        {input.fq} \
        -o {output.fa}"

# ------------------------------------------------------------------------------
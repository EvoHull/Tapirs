# ==================================================
# QUALITY CONTROL SNAKEFILE
# ==================================================
# --------------------------------------------------
# fastp, control for sequence quality and pair reads
# --------------------------------------------------

configfile: "config.yaml"

# ruleorder: vsearch_dereplication > empty_fasta_workaround > vsearch_rereplication

rule fastp_trim_and_merge:
    message:
        "Beginning fastp quality control of raw data"
    conda:
        "../envs/environment.yaml"
    input:
        read1 = "data/01_demultiplexed/{library}/{sample}.R1.fastq.gz",
        read2 = "data/01_demultiplexed/{library}/{sample}.R2.fastq.gz"
    output:
        out1 = "results/02_trimmed/{library}/{sample}.R1.fastq.gz",
        out2 = "results/02_trimmed/{library}/{sample}.R2.fastq.gz",
        out_unpaired1 = "results/02_trimmed/{library}/{sample}.unpaired.R1.fastq.gz",
        out_unpaired2 = "results/02_trimmed/{library}/{sample}.unpaired.R2.fastq.gz",
        out_failed = "results/02_trimmed/{library}/{sample}.failed.fastq.gz",
        merged = "results/02_trimmed/{library}/{sample}_merged.fastq.gz",
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
        --qualified_quality_phred {config[FASTP_qual_phred]} \
        --length_required {config[FASTP_len_required]} \
        --cut_tail \
        --trim_front1 {config[FASTP_trim_front1]} \
        --trim_front2 {config[FASTP_trim_front2]} \
        --max_len1 {config[FASTP_max_len1]} \
        --max_len2 {config[FASTP_max_len2]} \
        --merge \
        --merged_out {output.merged} \
        --overlap_len_require {config[FASTP_overlap_len]} \
        --correction \
        "


rule keep_fwd_unpaired:  # needs work
    input:
        merged = "results/02_trimmed/{library}/{sample}_merged.fastq.gz",
        out_unpaired1 = "results/02_trimmed/{library}/{sample}.unpaired.R1.fastq.gz"
    output:
<<<<<<< HEAD
        "results/02_trimmed/{library}/{sample}_merged.fastq.gz"
=======
        "results/02_trimmed/{library}/{sample}_catted.fastq.gz"
>>>>>>> f7332e27cf01cdd54302a7b421eabb5acd8014ca
    shell:
        "cat {input.out_unpaired1} {input.merged} > {output}"

# -----------------------------------------------------
# convert files from fastq to fasta
# -----------------------------------------------------

rule fastq_to_fasta:
    conda:
        "../envs/environment.yaml"
    input:
<<<<<<< HEAD
        "results/02_trimmed/{library}/{sample}_merged.fastq.gz"
    output:
        "results/02_trimmed/{library}/{sample}_merged.fasta",
=======
        "results/02_trimmed/{library}/{sample}_catted.fastq.gz"
    output:
        "results/02_trimmed/{library}/{sample}_catted.fasta",
>>>>>>> f7332e27cf01cdd54302a7b421eabb5acd8014ca
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
        "../envs/environment.yaml"
    input:
<<<<<<< HEAD
        "results/02_trimmed/{library}/{sample}_merged.fastq.gz"
=======
        "results/02_trimmed/{library}/{sample}_catted.fastq.gz"
>>>>>>> f7332e27cf01cdd54302a7b421eabb5acd8014ca
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
        "../envs/environment.yaml"
    input:
<<<<<<< HEAD
        "results/02_trimmed/{library}/{sample}_merged.fasta"
=======
        "results/02_trimmed/{library}/{sample}_catted.fasta"
>>>>>>> f7332e27cf01cdd54302a7b421eabb5acd8014ca
    output:
        "results/02_trimmed/{library}/{sample}_derep.fasta"
    shell:
        "vsearch \
        --derep_fulllength {input} \
        --sizeout \
        --minuniquesize {config[VSEARCH_minuniqsize]} \
        --output {output} \
        "


# rule empty_fasta_workaround:
#     input:
#         "results/02_trimmed/{library}/{sample}.merged.tmp.derep.fasta"
#     output:
#         denoise = "results/02_trimmed/{library}/{sample}.merged.derep.fasta",
#         rerep = "results/rereplicated/{library}/{sample}.fasta"
#     priority:
#         1
#     shell:
#         """
#         if wc -l {input} > 0
#         then
#             cp {input} {output.denoise}
#         else
#             cp {input} {output.rerep}
#         fi
#         """

# with open({input},'r') as file:
#     count=0
#     for line in file.readlines():
#         count+=1
#         if count > 0:
#             df=pd.read_csv(file.name, sep='\t', header=None,skiprows=1)
#             df[0]=pd.to_numeric(df[0].str.split(r"size=", expand=True)[1])
#         else:
#             df=pd.DataFrame([([0])+(['unidentified']*10)])


# -----------------------------------------------------
# denoise
# -----------------------------------------------------

rule vsearch_denoising:
    conda:
        "../envs/environment.yaml"
    input:
        "results/02_trimmed/{library}/{sample}_derep.fasta"
    output:
        fasta = "results/03_denoised/{library}/{sample}_denoise.fasta"
    #params:
    #    log="reports/denoise/{library}/vsearch.log"
    shell:
        """
        set +e
        vsearch --sizein --sizeout --cluster_unoise {input} --centroids {output.fasta} --minsize {config[VSEARCH_minsize]} --unoise_alpha {config[VSEARCH_unoise_alpha]}
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 0
        else
            exit 0
        fi
        """


# -----------------------------------------------------
# chimera removal
# -----------------------------------------------------

rule vsearch_dechimerisation: # output needs fixing
    conda:
        "../envs/environment.yaml"
    input:
        "results/03_denoised/{library}/{sample}_denoise.fasta"
    output: # fix
        text = "results/03_denoised/{library}/{sample}_chimera.txt",
        fasta = "results/03_denoised/{library}/{sample}_nc.fasta"
    params:
        db = config["dechim_blast_db"]
    shell:
        """
        set +e
        vsearch --uchime_ref {input} --db {params.db} --mindiffs {config[VSEARCH_mindiffs]} --mindiv {config[VSEARCH_mindiv]} --uchimeout {output.text} --nonchimeras {output.fasta}
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 0
        else
            exit 0
        fi
        """


# ------------------------------------------------------
# re-replication
# -------------------------------------------------------

rule vsearch_rereplication:
    conda:
        "../envs/environment.yaml"
    input:
        "results/03_denoised/{library}/{sample}_nc.fasta",
        # rule("empty_fasta_workaround")
    output:
        "results/rereplicated/{library}/{sample}_rerep.fasta"
    threads:
        6
    shell:
        """
        set +e
        vsearch --rereplicate {input} --output {output}
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 0
        else
            exit 0
        fi
        """

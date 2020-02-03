#-----------------------------------------------------
# metaclese-qc
# ---------
# A metabarcoding workflow for quality control
# before taxonomic assignment
#-----------------------------------------------------

# get the sequence files into a list
SAMPLES, = glob_wildcards("{sample}.R[1,2].fastq.gz")

#-----------------------------------------------------
rule all_qc:
    input:
        expand("data/{sample}.R[1,2].fastq.gz", sample=SAMPLES),
        expand("results/seqkit/{sample}.extendedFrags.fas", sample=SAMPLES),
        expand("reports/{sample}.seqkit_fastastats.md", sample=SAMPLES)

#-----------------------------------------------------
# demultiplex to gather data correctly
# this is a placeholder, needs fixing
#-----------------------------------------------------
rule demultiplex:
    input:
        raw_data="data/raw/{sample}.R[1,2].fastq.gz", # fix to take libraries
        dir_tsv="data/demultiplex/" #fix to take tsv files
    output:
        directory("results/demultiplex/")
    script:
        "scripts/demultiplex.py" # this script needs to be snakemakified

#-----------------------------------------------------
# fastp- this replaces trimmomatic etc
#-----------------------------------------------------
rule fastp:
    input:
        read1 = "data/{sample}.R1.fastq.gz",
        read2 = "data/{sample}.R2.fastq.gz"
    output:
        R1samples = "reports/fastpR1samples.txt",
        out1 = "results/fastp/{input.read1}",
        out2 = "results/fastp/{input.read2}"
    shell:
        "fastp --in1 {input.read1} --out1 {output.out1} --in2 {input.read2} --out2 {output.out2} | echo {input.read1} > {output.R1samples}.txt"

#-----------------------------------------------------
# denoise (replaces clustering)
#-----------------------------------------------------
rule denoise:
    input:
        "data/Lib1-Oct-SH1a.R1.fastq.gz"
    output:
        "denoised-sequences"
    shell:
        "vsearch --cluster_unoise {input} --alnout {output} --id 0.8 --iddef 0 --notrunclabels"

#-----------------------------------------------------
# flash, paired-end read merger
#-----------------------------------------------------
rule flash:
    input:
        forward="results/fastp/data/{sample}.R1.fastq.gz.out1.fastq",
        reverse="results/fastp/data/{sample}.R2.fastq.gz.out2.fastq"
        # does this find the same sample name or just A sample name?
        # I could get a list here then iterate for R1/2 ??
    params:
        prefix="{sample}",
        maxoverlap="150" # default is 65bp
    output:
        destination=directory("results/flash/{sample}"),
        pholder="results/flash/{sample}/extendedFrags.fastq.gz"
        #logfile="reports/{sample}_logfile.log" # this causes problems, why?
    shell:
        "flash -z -o {params.prefix} -M {params.maxoverlap} -d {output.destination} {input.forward} {input.reverse}"
        # | tee {output.logfile}" # I should try to write this to /reports?

#-----------------------------------------------------
# seqkit to convert files from fastq to fasta
#-----------------------------------------------------
rule seqkit:
        input:
            "results/flash/{sample}.extendedFrags.fastq.gz"
        output:
            directory("results/seqkit/{sample}"),
            "results/seqkit/{sample}/extendedFrags.fas"
        shell:
            "seqkit fq2fa {input} -o {output}"

#-----------------------------------------------------
# seqkit to write simple report on fasta files
#-----------------------------------------------------
rule seqkitstats:
        input:
            "results/seqkit/{sample}.fasta.gz"
        output:
            "reports/{sample}.seqkit_fastastats.md"
        shell:
            "seqkit stats {input} | csvtk csv2md -t -o {output}"

#-----------------------------------------------------
# tidy fastp report
#-----------------------------------------------------
rule fastp_report_mover: # moving reports from fastp
    input:
        "fastp.*" # check this syntax
    output:
        directory("reports/")
    shell:
        "mv {input} {output}"

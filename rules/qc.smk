#-----------------------------------------------------
# Tapirs-qc
# ---------
# A metabarcoding workflow for quality control
# before taxonomic assignment
#-----------------------------------------------------

# get the sequence files into a list
SAMPLES, = glob_wildcards("{sample}.R[1,2].fastq.gz")

#-----------------------------------------------------
rule all_qc:
    input:
        expand("data/demultiplexed/{sample}.R[1,2].fastq.gz", sample=SAMPLES),
        expand("results/seqkit/{sample}.extendedFrags.fas", sample=SAMPLES),
        expand("reports/{sample}.seqkit_fastastats.md", sample=SAMPLES)

#-----------------------------------------------------
# fastp- this replaces trimmomatic etc
#-----------------------------------------------------
rule fastp:
    input:
        read1 = "data/demultiplexed/{sample}.R1.fastq.gz",
        read2 = "data/demultiplexed/{sample}.R2.fastq.gz"
    output:
        R1samples = "reports/fastpR1samples.txt",
        out1 = "results/fastp/{input.read1}",
        out2 = "results/fastp/{input.read2}"
    shell:
        "fastp --in1 {input.read1} --out1 {output.out1} --in2 {input.read2} --out2 {output.out2} | echo {input.read1} > {output.R1samples}.txt"

#-----------------------------------------------------
# denoise (replaces clustering)
#       mostly pseudocode
#-----------------------------------------------------
rule denoise:
    input:
        "data/{sample}.fastq.gz"
    output:
        "denoised-sequences"
    shell:
        "vsearch --cluster_unoise {input} --alnout {output} --id 0.8 --iddef 0 --notrunclabels"

#-----------------------------------------------------
# flash, paired-end read merger
#-----------------------------------------------------
rule flash:
    input:
        forward="results/fastp/{sample}.R1.fastq.gz.out1.fastq",
        reverse="results/fastp/{sample}.R2.fastq.gz.out2.fastq"
    params:
        prefix="{sample}",
        maxoverlap="150" # default is 65bp
    output:
        destination=directory("results/flash/{sample}"),
        pholder="results/flash/{sample}/extendedFrags.fastq.gz"
    shell:
        "flash -z -o {params.prefix} -M {params.maxoverlap} -d {output.destination} {input.forward} {input.reverse}"
        # | tee {output.logfile}" # try to write this to /reports?

#-----------------------------------------------------
# seqkit to convert files from fastq to fasta
#-----------------------------------------------------
rule seqkit_fq2fa:
        input:
            "results/flash/{sample}.extendedFrags.fastq.gz"
        output:
            directory("results/seqkit/{sample}"),
            "results/seqkit/{sample}/extendedFrags.fas"
        shell:
            "seqkit fq2fa {input} -o {output}"

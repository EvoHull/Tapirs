#-----------------------------------------------------
# Tapirs-kraken2
# ---------
# A metabarcoding workflow to assign taxonomy using Kraken2
#-----------------------------------------------------

# get the sequence files into a list
SAMPLES, = glob_wildcards("{sample}.R[1,2].fastq.gz")

#-----------------------------------------------------
rule all:
    input:
        expand("data/{sample}.R[1,2].fastq.gz", sample=SAMPLES)

#-----------------------------------------------------
# Kraken2 assign taxonomy
#-----------------------------------------------------
rule kraken:
    input:
        query: "data/kraken/query/R2.fasta"
        database: directory("data/kraken/db/NCBI_nt")
      output:
        "results/kraken/R2out.txt"
    shell:
      "kraken2 --db {input.database} {input.query} --use-names --report {output}"

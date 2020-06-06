# ==================================================
# KRAKEN ANALYSIS
# ==================================================

configfile: "config.yaml"

# --------------------------------------------------
# Kraken, kmer based taxonomic id
# --------------------------------------------------

rule kraken2:
    conda:
        "../envs/environment.yaml"
    input:
        "results/07_rereplicated/{library}/{sample}.rerep.fasta"
    output:
        kraken_outputs = "results/kraken/outputs/{library}/{sample}.tsv",
        kraken_reports = "results/kraken/reports/{library}/{sample}.txt"
    threads:
        6
    params:
        confidence = "0.0",
        kraken_db = directory(config["kraken_db"]) # This is specified but not called in the shell command - Mike
    shell:
        "kraken2 \
        --db {params.kraken_db} {input} \
        --memory-mapping \
        --threads {threads} \
        --confidence {params.confidence} \
        --output {output.kraken_outputs} \
        --report {output.kraken_reports} \
        "

#-----------------------------------------------------
# Kraken output to BIOM format
#-----------------------------------------------------

rule kraken_biom_and_tsv:
    conda:
        "../envs/environment.yaml"
    input:
        expand("results/kraken/reports/{library}/{sample}.txt")
    output:
        "results/kraken/{my_experiment}.biom" #my_experiment=config["my_experiment"])
    params:
        "results/kraken/reports/*/*.txt"
    shell:
        "kraken-biom --max F -o {output} --fmt hdf5 {params}"
        "kraken-biom --max F -o {output} --fmt tsv {params}"

# #---------------------------------------------------
# # Biom convert, BIOM to tsv
# #---------------------------------------------------

# rule kraken_biom_to_tsv:
#     conda:
#         "../envs/environment.yaml"
#     input:
#         expand("results/kraken/{my_experiment}.biom", my_experiment=config["my_experiment"])
#     output:
#         expand("results/kraken/{my_experiment}.tsv", my_experiment=config["my_experiment"])
#     threads:
#         6
#     shell:
#         "biom convert -i {input} -o {output} --to-tsv --header-key taxonomy"

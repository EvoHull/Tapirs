# ==================================================
# KRAKEN2 ANALYSIS
# ==================================================

# configfile: "config.yaml"

# --------------------------------------------------
# Kraken, kmer based taxonomic id
# --------------------------------------------------

rule kraken2:
    conda:
        config['conda']
    input:
        reads = "results/09_rereplicated/{LIBRARIES}/{SAMPLES}.rerep.fasta"
    output:
        reports = "results/kraken2/reports/{LIBRARIES}/{SAMPLES}.txt",
        outputs = "results/kraken2/outputs/{LIBRARIES}/{SAMPLES}.krk"
    shell:
        "kraken2 --db {config[kraken2_db]} {input.reads} \
        --threads {config[kraken2_threads]} \
        --confidence {config[kraken2_confidence]} \
        --report {output.reports} \
        --output {output.outputs}"

#-----------------------------------------------------
# Kraken output to BIOM format
#-----------------------------------------------------

# rule kraken_biom_and_tsv:
#     conda:
#         "../envs/environment.yaml"
#     input:
#         expand("results/kraken/reports/{library}/{sample}.txt", sample=SAMPLES, library=LIBRARIES)
#     output:
#         biom = "results/kraken/{library}/{sample}.biom",
#         text = "results/kraken/reports/{library}/{sample}.txt"
#     shell:
#         """
#         kraken-biom --max F --fmt hdf5 -o {output.biom};
#         kraken-biom --max F --fmt tsv -o {output.text}
#         """

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

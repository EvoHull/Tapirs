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
        seqs = expand("results/07_rereplicated/{library}/{sample}.rerep.fasta", sample=SAMPLES, library=LIBRARIES),
        kraken_db = config["kraken_db"] # This is specified but not called in the shell command - Mike
    output:
        kraken_outputs = "results/kraken/outputs/{library}/{sample}.tsv",
        kraken_reports = "results/kraken/reports/{library}/{sample}.txt"
    threads:
        6
    params:
        confidence = "0.0"
    shell:
        "kraken2 \
        --db {input.kraken_db} {input.seqs} \
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
        expand("results/kraken/reports/{library}/{sample}.txt", sample=SAMPLES, library=LIBRARIES)
    output:
        biom = "results/kraken/{my_experiment}.biom",
        txt = "results/kraken/reports/*/*.txt"
    shell:
        "kraken-biom --max F -o {output.biom} --fmt hdf5 {output.txt}"
        "kraken-biom --max F -o {output.biom} --fmt tsv {output.txt}"

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


#-----------------------------------------------------
# Kraken taxonomic output to html
# produces interactive displays of taxonomic diversity
#-----------------------------------------------------
rule kraken_recentrifuge_fig
    conda:
        "../envs/environment.yaml"
    input:
        taxdb = config["taxdump"],
        krakenout = "results/kraken/outputs/{library}/{sample}.tsv"
    output:
        "reports/recentrifuge/{library}/{sample}.tsv.html"
    shell:
        "rcf -n {input.taxdb} -k {input.krakenout} -o {output}"
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
        reads = "results/09_rereplicated/{LIBRARIES}/{SAMPLES}.rerep.fasta"
    output:
        reports = "results/kraken/{LIBRARIES}/{SAMPLES}.txt",
        outputs = "results/kraken/{LIBRARIES}/{SAMPLES}.krk"

    threads:
        10
    shell:
        "kraken2 --db {config[kraken_db]} {input.reads} \
        --use-names --memory-mapping \
        --threads {threads} \
        --confidence 0 \
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


rule get_centrifuge:
    output:
        "recentrifuge/retaxdump"
    shell:
        """
        git clone https://github.com/khyox/recentrifuge.git
        """


#-----------------------------------------------------
# Recentrifuge produces interactive html displays
# of kraken output
#-----------------------------------------------------

rule kraken_recentrifuge_fig:
    conda:
        "../envs/environment.yaml"
    input:
        taxdb = config["taxdump"].replace("rankedlineage.dmp", ""),
        kraken_out_N = "results/kraken/{LIBRARIES}/{SAMPLES}.krk",
        rcf = "recentrifuge/retaxdump"
    output:
        "reports/recentrifuge/{LIBRARIES}/{SAMPLES}.krk.html"
    shell:
        "./recentrifuge/rcf -n {input.taxdb} -k {input.kraken_out_N} -o {output}"

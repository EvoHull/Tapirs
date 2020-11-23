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
        kraken_output = "results/kraken/outputs/{library}/{sample}.krk",
        kraken_report = "results/kraken/reports/{library}/{sample}.txt"
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
        --output {output.kraken_output} \
        --report {output.kraken_report} \
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
        biom = "results/kraken/{library}/{sample}.biom",
        text = "results/kraken/reports/{library}/{sample}.txt"
    shell:
        """
        kraken-biom --max F --fmt hdf5 -o {output.biom};
        kraken-biom --max F --fmt tsv -o {output.text}
        """

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
# Recentrifuge produces interactive html displays 
# of kraken output
#-----------------------------------------------------

rule kraken_recentrifuge_fig:
    conda:
        "../envs/environment.yaml"
    input:
        taxdb = config["taxdump"],
        # makes 1 report per library containing all samples, needs .krk extension
        # kraken_out = expand("results/kraken/outputs/{library}", library=LIBRARIES),
        kraken_out_N = expand("results/kraken/outputs/{library}/{sample}.krk", sample=SAMPLES, library=LIBRARIES)
    output:
        "reports/recentrifuge/{library}/{sample}.html"
    shell:
        "rcf -n {input.taxdb} -k {input.kraken_out_N} -o {output}"
# ==================================================
# SINTAX ANALYSIS
# ==================================================

configfile: "config.yaml"

# --------------------------------------------------
# sintax, kmer similarity taxonomic ID
# --------------------------------------------------

rule sintax:
    conda:
        "../envs/environment.yaml"
    input:
        query = "results/07_rereplicated/{library}/{sample}.rerep.fasta",
        database = config["sintax_db"]
    output:
        tsv = "results/sintax/{library}/{sample}.sintax.tsv",
    shell:
        "vsearch --sintax \
        {input.query} \
        --db {input.database} \
        --tabbedout {output.tsv} \
        --strand both \
        --sintax_cutoff {config[SINTAX_cutoff]} \
        "
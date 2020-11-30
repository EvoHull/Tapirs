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
        query = "results/09_rereplicated/{LIBRARIES}/{SAMPLES}.rerep.fasta",
        database = config["sintax_db"]
    output:
        tsv = "results/sintax/{LIBRARIES}/{SAMPLES}.sintax.tsv",
    shell:
        "vsearch --sintax \
        {input.query} \
        --db {input.database} \
        --tabbedout {output.tsv} \
        --strand both \
        --sintax_cutoff {config[SINTAX_cutoff]} \
        "

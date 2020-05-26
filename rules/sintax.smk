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
        query = "results/07_rereplicated/{library}/{sample}.rerep.fasta"
    output:
        "results/sintax/{library}/{sample}.reads.sintax"
    params:
        database = config["sintax_db"]
    shell:
        "vsearch -sintax \
        {input.query} \
        -db {params.database} \
        -tabbedout {output} \
        -strand both \
        -sintax_cutoff {config[SINTAX_cutoff]} \
        "

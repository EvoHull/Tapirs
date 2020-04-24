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
        query = "results/rereplicated/{library}/{sample}_rerep.fasta"
    output:
        "results/sintax/{library}/{sample}_reads.sintax"
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

# ==================================================
# BIOM
# ==================================================

configfile: "config.yaml"

# --------------------------------------------------
# biom
# --------------------------------------------------

rule biom:
    conda:
        "../envs/environment.yaml"
    input:
        "results/kraken.reports/{library}/{sample}.txt"
    output:

    threads:
        6
    shell:
        "kraken-biom \
        {input} \
        -o {output}"

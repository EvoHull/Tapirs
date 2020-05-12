# ==================================================
# BIOM
# ==================================================

configfile: "config.yaml"

# --------------------------------------------------
# Kraken, create biom format file from report
# all libraries should be pooled to 1 biom file
# --------------------------------------------------

rule kraken-to-biom:
    conda:
        "../envs/environment.yaml"
    input:
        "results/kraken/reports/{library}/{sample}.txt"
    output:
        "results/kraken/reports/kraken_biom.hdf5"
    threads:
        6
    shell:
        "kraken-biom {input} -o {output}"

# ---------------------------------------------------------
# MLCA, convert mlca tsv output to biom format
# ---------------------------------------------------------

rule mlca-to-biom:
    conda:
        "../envs/environment.yaml"
    input:
        "results/mlca/{library}/{sample}_lca.tsv"
    output:
        "results/mlca/mlca_biom.hdf5"
    shell:
        "biom convert -i {input} -o {output} --to-hdf5"


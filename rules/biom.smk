# ==================================================
# BIOM, file manipulations involving the BIOM format
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
        expand("results/kraken/reports/{sample.library}/{sample.sample}.txt",
               sample=sample.reset_index().itertuples()),
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
        expand("results/mlca/{sample.library}/{sample.sample}.lca.tsv"
               sample=sample.reset_index().itertuples()),
    output:
        "results/mlca/mlca_biom.hdf5"
    shell:
        "biom convert -i {input} -o {output} --to-hdf5"

# ---------------------------------------------------------
# sintax, convert tabbedout to BIOM hdf5
# ---------------------------------------------------------

rule sintax_tsv_to_BIOM:
    input:
        expand("reports/sintax/{sample.library}/{sample.sample}.sintax.tsv"
               sample=sample.reset_index().itertuples()),
    output:
        "reports/sintax/{library}/{sample}.sintax.biom"
    shell:
        """
        biom convert {input} {output} --table-type="OTU table" --to-hdf5
        """
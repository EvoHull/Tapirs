# ==================================================
# BIOM, file manipulations involving the BIOM format
# ==================================================

configfile: "config.yaml"

# ---------------------------------------------------------
# MLCA, convert mlca tsv output to biom format
# ---------------------------------------------------------

rule mlca_to_biom:
    conda:
        "../envs/environment.yaml"
    input:
        "results/mlca/{LIBRARIES}/{SAMPLES}.lca.tsv"
    output:
        "results/mlca/mlca_biom.hdf5"
    shell:
        "biom convert -i {input} -o {output} --to-hdf5"

# ---------------------------------------------------------
# sintax, convert tabbedout to BIOM hdf5
# ---------------------------------------------------------

rule sintax_tsv_to_BIOM:
    conda:
        "../envs/environment.yaml"
    input:
        "reports/sintax/{LIBRARIES}/{SAMPLES}.sintax.tsv"
    output:
        "reports/sintax/{LIBRARIES}/{SAMPLES}.sintax.biom"
    shell:
        """
        biom convert {input} {output} --table-type="OTU table" --to-hdf5
        """

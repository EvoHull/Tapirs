# ==================================================
# PCA Plot
# ==================================================

configfile: "config.yaml"

# --------------------------------------------------
# Convert Kraken2 biom files to vegan format
# --------------------------------------------------

rule vegan_format_conversion:
    conda:
        "../envs/environment.yaml"
    input:
        expand("results/kraken/{my_experiment}.tsv", my_experiment=config["my_experiment"])
    output:
        expand("results/kraken/{my_experiment}_vegan.tsv", my_experiment=config["my_experiment"])
    shell:
        "sed -i "/OTU/s/^\#//g" {output}"
        # Need to run this to include sample names in the tsv table

# --------------------------------------------------
# Produce PCA plot
# --------------------------------------------------

rule pca_plot:
    conda:
        "../envs/environment.yaml"
    input:
        kraken = expand("results/kraken/{my_experiment}.tsv", my_experiment=config["my_experiment"])
        blast = "results/" + config['my_experiment'] + "blast" + config['MLCA_identity'] + ".tsv"
    output:
        "reports/r_analyses/pca_plot.png"
    script:
        "scripts/vegan_analyses_test.R"

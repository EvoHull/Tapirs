# ==============================================================================
# TAPIRS
# A reproducible metabarcoding workflow using snakemake
# ==============================================================================

import pandas as pd
configfile: "config.yaml"

# --------------------------------------------------
# Wildcarding library and sample
# --------------------------------------------------

library = pd.read_csv(config["libraries"], sep='\t',
                      header=0, index_col="Library")
# library = pd.read_table(config["libraries"], index_col="library")
sample = pd.read_csv(config["samples"], sep='\t',
                     header=0, index_col=0, dtype=str)
# sample = pd.read_table(config["samples"], index_col=[
#                        "Library", "Sample"], dtype=str)
sample.index = sample.index.set_levels(
    [i.astype(str) for i in sample.index.levels])

# ---------------------------------------------------
# Target rules
# ---------------------------------------------------
rule all:
    input:
# results
        expand("results/kraken/outputs/{sample.library}/{sample.sample}.tsv",
               sample=sample.reset_index().itertuples()),
        expand("results/kraken/reports/{sample.library}/{sample.sample}.txt",
               sample=sample.reset_index().itertuples()),
# snakemake reports
        "reports/rulegraph_dag.svg",
        "reports/rulegraph_dag.png",
        "reports/snakemake-report.html",
        expand("reports/archived_envs/{conda_envs}",
               conda_envs=config["conda_envs"]),
# fastp and multiqc reports
        expand("reports/fastp/{sample.library}/{sample.sample}_fastp.json",
               sample=sample.reset_index().itertuples()),
        expand("reports/fastp/{sample.library}/{sample.sample}_fastp.html",
               sample=sample.reset_index().itertuples()),
        expand("reports/multiqc/{library}.multiqc.html"),
# vsearch reports
        #expand("reports/vsearch/{sample.library}/{sample.sample}.denoise.biom", sample=sample.reset_index().itertuples()),
        expand("reports/vsearch/{sample.library}/{sample.sample}_fq_eestats",
               sample=sample.reset_index().itertuples()),
        expand("reports/vsearch/{sample.library}/{sample.sample}_fq_readstats",
               sample=sample.reset_index().itertuples()),
# kraken and sintax reports
        #expand("results/kraken/{my_experiment}.biom", my_experiment=config["my_experiment"]),
        expand("results/kraken/{my_experiment}.tsv",
               my_experiment=config["my_experiment"]),
        expand("results/sintax/{sample.library}/{sample.sample}_reads.sintax",
               sample=sample.reset_index().itertuples()),
# ????
        expand("reports/{my_experiment}.tsv",
               my_experiment=config["my_experiment"]),

# -----------------------------------------------------
# Rule files
# -----------------------------------------------------
include: "rules/qc.smk"
include: "rules/blast.smk"
include: "rules/kraken.smk"
include: "rules/mlca.smk"
include: "rules/sintax.smk"
include: "rules/reports.smk"

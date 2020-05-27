# ==============================================================================
# TAPIRS
# A reproducible metabarcoding workflow using snakemake
# ==============================================================================

import pandas as pd
configfile: "config.yaml"

# --------------------------------------------------
# Wildcarding library and sample
# --------------------------------------------------

library = pd.read_table(config["libraries"], index_col="library")
# sample = pd.read_table(config["samples"], dtype=str)
sample = pd.read_table(config["samples"], index_col=["library", "sample"], dtype=str)
sample.index = sample.index.set_levels([i.astype(str) for i in sample.index.levels])


# library = pd.read_csv(config["libraries"], sep='\t',
#                      header=0, index_col="library")
# # library = pd.read_table(config["libraries"], index_col="library")
# sample = pd.read_csv(config["samples"], sep='\t',
#                      header=0, index_col=["library", "sample"], dtype=str)
# # sample = pd.read_table(config["samples"], index_col=[
# #                        "library", "sample"], dtype=str)
# sample.MultiIndex = sample.MultiIndex.set_levels(
#     [sample.i.astype(str) for i in sample.MultiIndex.levels])

# sample.index = sample.index.set_levels(
#     [i.astype(str) for i in sample.index.levels])
# df.index = df.index.set_levels(df.index.levels[2].astype(int), level=2)
#  so `i.astype(str)` should be `sample.i.astype(str)`


# ---------------------------------------------------
# Target rules
# ---------------------------------------------------
rule all:
    input:
# snakemake reports and envs
       #  "report.html",
        "reports/rulegraph_dag.png",
       #  "snakemake-report.html",
        expand("reports/archived_envs/{conda_envs}",
               conda_envs=config["conda_envs"]),
# fastp, multiQC, and seqkit reports
        expand("reports/fastp/{sample.library}/{sample.sample}.fastp.json",
               sample=sample.reset_index().itertuples()),
        expand("reports/fastp/{sample.library}/{sample.sample}.fastp.html",
               sample=sample.reset_index().itertuples()),
        # expand("reports/multiqc/{library}.multiqc.html", 
        #        library=library.reset_index().itertuples()),
        expand("reports/seqkit/{sample.library}.trimmed.seqkit-stats.tsv", 
               sample=sample.reset_index().itertuples()),
        expand("reports/seqkit/{sample.library}.trimmed.seqkit-stats.md", 
               sample=sample.reset_index().itertuples()),
        expand("reports/seqkit/{sample.library}.concat.seqkit-stats.tsv",
               sample=sample.reset_index().itertuples()),
        expand("reports/seqkit/{sample.library}.concat.seqkit-stats.md",
               sample=sample.reset_index().itertuples()),  
       #  expand("results/03_merged/{sample.library}/{sample}.concat.fasta",
       #         sample=sample.reset_index().itertuples()),     
# vsearch reports
        #expand("reports/vsearch/{sample.library}/{sample.sample}.denoise.biom", 
       #       sample=sample.reset_index().itertuples()),
        expand("reports/vsearch/{sample.library}/{sample.sample}.concat.fq_eestats",
               sample=sample.reset_index().itertuples()),
        expand("reports/vsearch/{sample.library}/{sample.sample}.concat.fq_readstats",
               sample=sample.reset_index().itertuples()),
# sintax
        expand("results/sintax/{sample.library}/{sample.sample}.reads.sintax",
               sample=sample.reset_index().itertuples()),
# blast mlca
        expand("results/blast/{sample.library}/{sample.sample}.blast.tsv",
               sample=sample.reset_index().itertuples()),  # optional
        expand("results/mlca/{sample.library}/{sample.sample}.lca.tsv",
               sample=sample.reset_index().itertuples()),
# mlca-tsv.py
       #  expand("reports/{my_experiment}.tsv",
       #         my_experiment=config["my_experiment"]),
# Kraken
        expand("results/kraken/outputs/{sample.library}/{sample.sample}.tsv",
               sample=sample.reset_index().itertuples()),
        expand("results/kraken/reports/{sample.library}/{sample.sample}.txt",
               sample=sample.reset_index().itertuples()),
        expand("results/kraken/{my_experiment}.tsv",
              my_experiment=config["my_experiment"]),
       #  expand("results/kraken/{my_experiment}.biom", 
       #        my_experiment=config["my_experiment"]),

# -----------------------------------------------------
# Rule files
# -----------------------------------------------------
include: "rules/qc.smk"
include: "rules/blast.smk"
include: "rules/kraken.smk"
include: "rules/lca.smk"
include: "rules/sintax.smk"
include: "rules/reports.smk"

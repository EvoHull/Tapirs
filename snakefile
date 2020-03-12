###########################################################################################################################

                                                ######   TAPIRS   #####
                                # A reproducible metabarcoding workflow using snakemake

############################################################################################################################
## Setup

# Config file
configfile: "config.yaml"

# Reporting
report: "reports/tapirs.rst" # Flag "$ snakemake" with "--report" to use

###########################################################################################################################
## Wildcarding library and sample

import pandas as pd

library = pd.read_table(config["libraries"], index_col="library")
sample = pd.read_table(config["samples"], index_col=["library", "sample"], dtype=str)
sample.index = sample.index.set_levels([i.astype(str) for i in sample.index.levels])

############################################################################################################################
#-------------------------------------------------------------------------------
# Target rules
rule all:
    input:
# krona taxonomy database
        #"data/databases/krona/taxonomy.tab",
# results ----------------------------------------------------------------------
        expand("results/kraken/outputs/{sample.library}/{sample.sample}.tsv", sample=sample.reset_index().itertuples()),
        expand("results/kraken/reports/{sample.library}/{sample.sample}.txt", sample=sample.reset_index().itertuples()),

# reports ----------------------------------------------------------------------
        expand("reports/fastp/{sample.library}/{sample.sample}.json", sample=sample.reset_index().itertuples()),
        expand("reports/fastp/{sample.library}/{sample.sample}.html", sample=sample.reset_index().itertuples()),
        #expand("reports/vsearch/{sample.library}/{sample.sample}.denoise.biom", sample=sample.reset_index().itertuples()),
        expand("reports/vsearch/{sample.library}/{sample.sample}_fq_eestats", sample=sample.reset_index().itertuples()),
        expand("reports/vsearch/{sample.library}/{sample.sample}_fq_readstats", sample=sample.reset_index().itertuples()),
        expand("reports/archived_envs/{conda_envs}", conda_envs=config["conda_envs"]),
        #expand("results/kraken/{my_experiment}.biom", my_experiment=config["my_experiment"]),
        expand("results/kraken/{my_experiment}.tsv", my_experiment=config["my_experiment"]),
        expand("reports/krona/kraken/{sample.library}/{sample.sample}.2.html", sample=sample.reset_index().itertuples()),
        expand("reports/krona/mlca/{sample.library}/{sample.sample}.html", sample=sample.reset_index().itertuples()),
        expand("results/sintax/{sample.library}/{sample.sample}_reads.sintax", sample=sample.reset_index().itertuples()),
        expand("reports/krona/sintax/{sample.library}/{sample.sample}.sintax.html", sample=sample.reset_index().itertuples()),
        expand("reports/mlca/mlca2tsv/{my_experiment}.tsv", my_experiment=config["my_experiment"]),
# for testing
        #expand("results/blast/{sample.library}/{sample.sample}_blast.taxed.out", sample=sample.reset_index().itertuples()),
        #expand("results/mlca/{sample.library}/{sample.sample}_lca.tsv", sample=sample.reset_index().itertuples()),
        #expand("results/02_trimmed/{sample.library}/{sample.sample}.{R}.fastq.gz", sample=sample.reset_index().itertuples(), R=config["R"]),
        #expand("results/02_trimmed/{sample.library}/{sample.sample}.unpaired.{R}.fastq.gz", sample=sample.reset_index().itertuples(), R=config["R"]),
        #expand("results/02_trimmed/{sample.library}/{sample.sample}.merged.fastq.gz", sample=sample.reset_index().itertuples()),
        #expand("results/03_denoised/{sample.library}/{sample.sample}.fasta", sample=sample.reset_index().itertuples(), R=config["R"]),
        #expand("results/blast/{sample.library}/{sample.sample}_blast.out", sample=sample.reset_index().itertuples()),
        #expand("reports/multiqc/{sample.library}.multiqc.html", sample=sample.reset_index().itertuples())
#-----------------------------------------------------
# Rule files
#-----------------------------------------------------
include: "rules/qc.smk"
include: "rules/blast.smk"
include: "rules/kraken.smk"
include: "rules/mlca.smk"
include: "rules/sintax.smk"
include: "rules/reports.smk"

##################################################################################################

# rule kt_taxonomy:
#     conda:
#         "envs/environment.yaml"
#     output:
#         "data/databases/krona/taxonomy.tab"
#     params:
#         "data/databases/krona/"
#     priority:
#         50
#     shell:
#         "rm -rf {params} \
#         && mkdir {params} \
#         && ktUpdateTaxonomy.sh {params}"

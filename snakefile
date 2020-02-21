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
# results ----------------------------------------------------------------------
        expand("results/02_trimmed/{library}/{sample}.{R}.fastq.gz", library=library, sample=sample, R=config["R"]),
        expand("results/02_trimmed/{library}/{sample}.unpaired.{R}.fastq.gz", library=library, sample=sample, R=config["R"]),
        expand("results/02_trimmed/{library}/{sample}.merged.fastq.gz", library=library, sample=sample),
        expand("results/03_denoised/{library}/{sample}.fasta", library=library, sample=sample, R=config["R"]),
        expand("results/blast/{library}/{sample}_blast.out", library=library, sample=sample),
        #expand("results/blast/{library}/{sample}_blast.taxed.out", library=library, sample=sample),
        expand("results/mlca/{library}/{sample}_lca.tsv", library=library, sample=sample),
# reports ----------------------------------------------------------------------
        expand("reports/fastp/{library}/{sample}.json", library=library, sample=sample),
        expand("reports/fastp/{library}/{sample}.html", library=library, sample=sample),
        expand("reports/vsearch/{library}/{sample}.denoise.biom", library=library, sample=sample),
        expand("reports/vsearch/{library}/{sample}_fq_eestats", library=library, sample=sample),
        expand("reports/vsearch/{library}/{sample}_fq_readstats", library=library, sample=sample),
        expand("reports/archived_envs/{conda_envs}", conda_envs=config["conda_envs"]),
        expand("results/kraken/{my_experiment}.tsv", my_experiment=config["my_experiment"]),
        expand("reports/krona/kraken/{sample.library}/{sample.sample}.html", sample=sample.reset_index().itertuples()),
        expand("reports/krona/mlca/{library}/{sample}.html", library=library, sample=sample)

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

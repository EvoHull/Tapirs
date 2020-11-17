# ==============================================================================
# TAPIRS
# A reproducible metabarcoding workflow using snakemake
# ==============================================================================

import pandas as pd
configfile: "config.yaml"

# --------------------------------------------------
# Load library and sample information
# --------------------------------------------------

samples_df = pd.read_table('samples.tsv').set_index("sample", drop=False)
SAMPLES = list(samples_df['sample'])
libraries_df = pd.read_table('libraries.tsv').set_index("library", drop=False)
LIBRARIES = list(libraries_df['library'])

# samples_df = pd.read_table(config['samples'], index_col="sample", header=0, drop=FALSE)
# SAMPLES = list(samples_df.index)
# libraries_df = pd.read_table(config['libraries'], index_col="library", header=0, drop=FALSE)
# LIBRARIES = list(libraries_df.index)

# ---------------------------------------------------
# Target rules
# ---------------------------------------------------

rule all:
    input:
        "reports/rulegraph_dag.png",
        expand("reports/archived_envs/{conda_envs}", 
                conda_envs=config["conda_envs"]),
# fastp, multiQC, and seqkit reports
        expand("reports/fastp/{library}/{sample}.fastp.json",
                sample=SAMPLES, library=LIBRARIES),
        expand("reports/fastp/{library}/{sample}.fastp.html",
                sample=SAMPLES, library=LIBRARIES),
        # expand("reports/seqkit/{library}/{sample}.trimmed.seqkit-stats.tsv", 
        #         sample=SAMPLES, library=LIBRARIES),
        # expand("reports/seqkit/{library}/{sample}.trimmed.seqkit-stats.md", 
        #         sample=SAMPLES, library=LIBRARIES),
        expand("reports/seqkit/{library}.concat.seqkit-stats.tsv",
                sample=SAMPLES, library=LIBRARIES),
        expand("reports/seqkit/{library}.concat.seqkit-stats.md",
                sample=SAMPLES, library=LIBRARIES),
        # expand("reports/multiqc/{library}.multiqc.html", 
        #        library=library.reset_index().itertuples()),
# vsearch reports
        # expand("reports/vsearch/{library}/{sample}.concat.fq_eestats",
        #         sample=SAMPLES, library=LIBRARIES),
        # expand("reports/vsearch/{library}/{sample}.concat.fq_readstats",
        #         sample=SAMPLES, library=LIBRARIES),
# sintax
        expand("results/sintax/{library}/{sample}.sintax.tsv",
                sample=SAMPLES, library=LIBRARIES),
# blast mlca
        expand("results/blast/{library}/{sample}.blast.tsv",
                sample=SAMPLES, library=LIBRARIES),
        expand("results/mlca/{library}/{sample}.lca.tsv",
                sample=SAMPLES, library=LIBRARIES),
# Kraken
        expand("results/kraken/outputs/{library}/{sample}.tsv",
                sample=SAMPLES, library=LIBRARIES),
        expand("results/kraken/reports/{library}/{sample}.txt",
                sample=SAMPLES, library=LIBRARIES),
        # expand("results/kraken/{my_experiment}.tsv", 
                # my_experiment=config["my_experiment"]),
        # expand("results/kraken/{my_experiment}.biom", 
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
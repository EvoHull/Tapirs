# ==============================================================================
# TAPIRS
# A reproducible metabarcoding workflow using snakemake
# ==============================================================================

import pandas as pd
configfile: "config.yaml"
report: "reports/snakemake-report.rst"
# --------------------------------------------------
# Load library and sample information
# --------------------------------------------------

with open('samples.tsv') as infile:
    sample_list = []
    for line in infile:
        sample_list.append(line.strip().replace('\t', '/'))
SAMPLES = sample_list

# ---------------------------------------------------
# Target rule
# ---------------------------------------------------

rule all:
    input:
# Final csv
        "results/"+config['my_experiment']+"blast"+config['MLCA_identity']+".tsv",
        expand("results/kraken/{sample}.txt", sample = SAMPLES),
# Reports
        "reports/dag_rulegraph.png",
# Archives
        "reports/archived_envs/tapirs.yaml",



# -----------------------------------------------------
# Rule files
# -----------------------------------------------------

include: "rules/qc.smk"
include: "rules/blast.smk"
include: "rules/kraken2.smk"
include: "rules/lca.smk"
include: "rules/vsearch.smk"
# include: "rules/sintax.smk"
include: "rules/reports.smk"

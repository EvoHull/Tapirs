# ==============================================================================
# TAPIRS
# A reproducible metabarcoding workflow using snakemake
# ==============================================================================

import pandas as pd

configfile: "config.yaml"
report: "reports/snakemake-report.rst"

# --------------------------------------------------
# LOAD LIBRARY AND SAMPLE INFORMATION


libraries_df = pd.read_table('output.tsv', header = None)  # Pull in libraries and samples from tsv

# Generate library wildcards
LIBRARIES = list(libraries_df[0])  # create list from column of dataframe
LIBRARIES = list(dict.fromkeys(LIBRARIES))  # dipping it into a dictionary to remove duplicates

# Generate samples wildcards
SAMPLES = list(libraries_df[1])  # create list from column of dataframe
SAMPLES = list(dict.fromkeys(SAMPLES))  # dipping it into a dictionary to remove duplicates

# Generate list of legitimate combinations
with open('output.tsv') as infile:
    real_combos = []
    for line in infile:
        real_combos.append(line.strip().replace('\t', '/'))


rule all:
    input:
        expand("results/sintax/{real_combos}.sintax.tsv", real_combos = real_combos),
        "results/kraken/{real_combos}.txt", real_combos = real_combos),
        "results/"+config["my_experiment"]+"blast"+config['MLCA_identity']+".tsv",

        "reports/dag_rulegraph.png",
        expand("reports/fastp/{real_combos}.fastp.html", real_combos = real_combos),
        expand("reports/recentrifuge/{real_combos}.krk.html", real_combos = real_combos),

        "reports/archived_envs/tapirs.yaml"

        # expand("results/mlca/{real_combos}.mlca_biom.hdf5", real_combos = real_combos),
        # expand("reports/sintax/{real_combos}.sintax.biom", real_combos = real_combos)


# -----------------------------------------------------
# RULE FILES

include: "rules/kraken2.smk"
include: "rules/sintax.smk"
include: "rules/qc.smk"
include: "rules/blast.smk"
include: "rules/lca.smk"
include: "rules/vsearch.smk"
include: "rules/reports.smk"
include: "rules/biom.smk"

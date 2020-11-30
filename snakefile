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

libraries_df = pd.read_table('output.tsv', header = None)

# Generate library wildcards
LIBRARIES = list(libraries_df[0])
LIBRARIES = list(dict.fromkeys(LIBRARIES))

# Generate samples wildcards
SAMPLES = list(libraries_df[1])
SAMPLES = list(dict.fromkeys(SAMPLES))

# Generate list of legitimate combinations
with open('output.tsv') as infile:
    real_combos = []
    for line in infile:
        real_combos.append(line.strip().replace('\t', '/'))


# ---------------------------------------------------
# Target rule
# ---------------------------------------------------

rule all:
    input:
# Sintax
        # expand("results/sintax/{real_combos}.sintax.tsv", real_combos = real_combos),
# Final csv
        "results/"+config['my_experiment']+"blast"+config['MLCA_identity']+".tsv",
        expand("results/kraken/{real_combos}.txt", real_combos = real_combos),
# Reports
        "reports/dag_rulegraph.png",
        expand("reports/fastp/{real_combos}.fastp.html", real_combos = real_combos),
# Archives
        "reports/archived_envs/tapirs.yaml",
# Biom
        expand("results/mlca/{real_combos}.mlca_biom.hdf5", real_combos = real_combos),
        # expand("reports/sintax/{real_combos}.sintax.biom", real_combos = real_combos)

# -----------------------------------------------------
# Rule files
# -----------------------------------------------------

include: "rules/qc.smk"
include: "rules/blast.smk"
include: "rules/kraken2.smk"
include: "rules/lca.smk"
include: "rules/vsearch.smk"
include: "rules/sintax.smk"
include: "rules/reports.smk"
include: "rules/biom.smk"

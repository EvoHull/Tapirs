# https://github.com/snakemake-workflows/rna-seq-star-deseq2/blob/master/Snakefile

import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

# we need a libraries.tsv and units.tsv files at the top level
# config.yaml then just lists them ie
units: units.tsv
libraries: libraries.tsv

# Note 'samples' changed to 'libraries'

libraries = pd.read_table(config["libraries"]).set_index("sample", drop=False)
validate(libraries, schema="schemas/libraries.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="schemas/units.schema.yaml")


##### target rules #####

rule all:
    input:
        expand(["results/diffexp/{contrast}.diffexp.tsv",
                "results/diffexp/{contrast}.ma-plot.svg"],
               contrast=config["diffexp"]["contrasts"]),
        "results/pca.svg",
        "qc/multiqc_report.html"

##### setup report #####

report: "report/workflow.rst"

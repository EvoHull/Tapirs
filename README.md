![tapirs_logo](documentation/docs/images/tapirs_seq.png)
# Tapirs
a snakemake workflow for reproducible metabarcoding

Tapirs is a reproducible modular workflow for the analysis of metabarcoding data.

Tapirs is compartmentalised into several modules, all contained in the rules/ directory. Each module performs a step of the workflow.

Each rule is assigned to a Conda environment containing the appropriate packages needed to perform. Using Conda ensures version control and prevents workflow failure through package incompatability.

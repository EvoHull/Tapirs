![tapirs_logo](docs/images/tapirs_seq.png)

# Tapirs

Tapirs is a reproducible modular workflow for the analysis of DNA metabarcoding data.

Tapirs uses the Snakemake workflow manager and is compartmentalised into several modules, each performing a step of the workflow. Tapirs is designed to be experimental, allowing you to test the effect of different approaches to data analysis.

Rules make use of Conda environments containing the appropriate packages needed to perform. Using Conda ensures version control and prevents workflow failure through package incompatability.

Tapirs was created by the EvoHull group, the University of Hull, UK

Instructions for installation, setup, and modification are contained within the [Tapirs documentation](https://tapirs.readthedocs.io)

## DAG overview of a workflow

![workflow graph](docs/images/dag.svg)

## Quickstart

Tapirs is curently in alpha/beta-release, not all features are present and not all bugs have been caught

Detailed instructions are given in the [Tapirs documentation](https://tapirs.readthedocs.io).

1. install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) (miniconda)
2. git clone the Tapirs repository
    - `git clone https://github.com/davelunt/Tapirs`
3. install snakemake in your base conda environment
    - `conda activate base`
    - `conda install -c bioconda -c conda-forge snakemake`
4. Place all library directories within the "data/01_demultiplexed/" directory ensuring they follow the format:
`data/01_demultiplexed/<library>/<sample>.<read>.fastq.gz`
5. Run the script to create the library and sample lists from your data.
`bash scripts/wildcarding.sh`
5. dry run `snakemake -s snakefile --use-conda --printshellcmds -n` to identify any issues
6. run `snakemake -s snakefile --use-conda --printshellcmds`

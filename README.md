# Tapirs Workflow

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/cc_tapirs.svg?branch=master)](https://travis-ci.org/snakemake-workflows/cc_tapirs)

Tapirs is a reproducible modular workflow for the analysis of DNA metabarcoding data.

Tapirs uses the Snakemake workflow manager and is compartmentalised into several modules, each performing a step of the workflow. Tapirs is designed to be experimental, allowing you to test the effect of different approaches to data analysis.

Rules make use of Conda environments containing the appropriate packages needed to perform. Using Conda ensures version control and prevents workflow failure through package incompatability.

Tapirs was created by the EvoHull group, the University of Hull, UK

Instructions for installation, setup, and modification are contained within the [Tapirs documentation](https://tapirs.readthedocs.io)

## DAG overview of a workflow

![workflow graph](docs/images/dag.svg)

## Authors

* Dave Lunt (@davelunt)
* Graham Sellers (@Graham-Sellers)
* Mike Winter (@mrmrwinter)
* Merideth Freiheit (@merfre)
* et al.

## Quickstart

Tapirs is curently in alpha/beta-release, not all features are present and not all bugs have been caught

Detailed instructions are given in the [Tapirs documentation](https://tapirs.readthedocs.io).

1. install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) (miniconda)
2. git clone the Tapirs repository
    - `git clone https://github.com/davelunt/Tapirs`
3. install snakemake
    - `conda install -c bioconda -c conda-forge snakemake`
4. Place all library directories within the "resources/libraries" directory ensuring they follow the format:
`libraries/01_demultiplexed/<library>/<sample>.<read>.fastq.gz`
5. Create the library and sample lists from your data.
5. Dry run `snakemake -npr` to identify any issues
6. Run `snakemake`

## Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `samples.tsv` to specify your sample setup.
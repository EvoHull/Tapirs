![tapirs_logo](documentation/docs/images/tapirs_seq.png)
# Tapirs

Tapirs is a reproducible modular workflow for the analysis of DNA metabarcoding data.

Tapirs uses the Snakemake workflow manager and is compartmentalised into several modules, all contained in the rules/ directory. Each module performs a step of the workflow.

Each rule is assigned to a Conda environment containing the appropriate packages needed to perform. Using Conda ensures version control and prevents workflow failure through package incompatability.

Tapirs was created by the EvoHull group, the University of Hull, UK

Instructions for installation, setup, and modification are contained within the [Tapirs documentation](documentation/docs/site/index.html)

DAG overview of a workflow
![workflow graph](documentation/docs/images/dag.svg)

## Instructions

#### CONDA and Snakemake

First, ensure a version of conda is intalled on your machine. This can be either anaconda or miniconda.
Miniconda is available at: https://docs.conda.io/en/latest/miniconda.html

Next, create and activate the tapirs conda environment.
This can be done with the following command:

$ conda create --file envs/tapirs.yaml

$ conda activate tapirs


### Database setup

#### BLAST
BLAST requires a local custom database for the workflow to run efficiently.
######## GS or DL - can you put a little bit in here about making blast databases etc - Mike ####


#### KRONA
Krona requires that its taxonomy database is updated/created
To do this, run the following commands, making sure you are within the Tapirs/ directory and the 'tapirs' conda environment.

$ mkdir data/databases/krona/

$ ktUpdateTaxonomy.sh data/databases/krona/


#### KRAKEN
Kraken needs its database to run.
######## GS, cna you sort this out please? - Mike   ########




### Input wrangling

Place all library directories within the "data/" directory (or edit path above), ensuring they follow the format "<library>.<sample>.<read>.fastq.gz".


### Wildcard generation

<insert command to run script to create tsvs for samples>

$ bash wildcarding.sh


### Config setup

Open the configuration file "config.yaml" and follow steps there.

Once this file has been configured we are ready to run the workflow.


#### Running the workflow

The workflow can be ran with the following line of code:

$ snakemake -s snakefile --use-conda --printshellcmds -n

The above line will start a dryrun. Assess the joblist of the dryrun and the fnal output files to ensure that your libraries and samples are being detected by the iteration steps.

If all looks well, run the workflow for real with the following code:

$ snakemake -s snakefile --use-conda --printshellcmds

![tapirs_logo](documentation/docs/images/tapirs_seq.png)
# Tapirs

Tapirs is a reproducible modular workflow for the analysis of DNA metabarcoding data.

Tapirs uses the Snakemake workflow manager and is compartmentalised into several modules, all contained in the rules/ directory. Each module performs a step of the workflow.

Each rule is assigned to a Conda environment containing the appropriate packages needed to perform. Using Conda ensures version control and prevents workflow failure through package incompatability.

Tapirs was created by the EvoHull group, the University of Hull, UK

Instructions for installation, setup, and modification are contained within the [Tapirs documentation](documentation/site/index.html)

## DAG overview of a workflow

![workflow graph](documentation/docs/images/dag.svg)

## Instructions

#### CONDA and Snakemake

First, ensure a version of conda is intalled on your machine. This can be either anaconda or miniconda.
Miniconda is available at: https://docs.conda.io/en/latest/miniconda.html

Once conda is installed, activate your base conda environment and install snakemake within it.
This can be done with the following command:

`conda activate base`

`conda install -c bioconda -c conda-forge snakemake`

Snakemake must be ran with the --use-conda flag.

Run the workflow from within this environment.

### Database setup

#### BLAST
BLAST requires a local custom database for the workflow to run efficiently.
######## GS or DL - can you put a little bit in here about making blast databases etc - Mike ####


#### KRONA
Kronas taxonomy database is updated/created from within the snakemake workflow.

#### KRAKEN
Kraken needs either the full Kraken database or a custom database to run.
######## GS, cna you sort this out please? - Mike   ########


#### SINTAX
Sintax also needs its own database.


### Input wrangling

The workflow assumes that input data has already been demultiplexed.

Place all library directories within the "data/01_demultiplexed/" directory (or edit path above), ensuring they follow the format,
`data/01_demultiplexed/<library>/<sample>.<read>.fastq.gz`


### Wildcard generation
Once input files have been placed in the `01_demultiplexed` directory, run the following script to create the library and sample lists.

`bash scripts/wildcarding.sh`


### Config setup

Open the configuration file `config.yaml` and follow steps there.

Once this file has been configured we are ready to run the workflow.


#### Running the workflow

The workflow can be ran with the following line of code:

`snakemake -s snakefile --use-conda --printshellcmds -n`

The above line will start a dryrun. Assess the joblist of the dryrun and the fnal output files to ensure that your libraries and samples are being detected by the iteration steps.

If all looks well, run the workflow for real with the following code:

`snakemake -s snakefile --use-conda --printshellcmds`

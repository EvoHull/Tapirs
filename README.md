![tapirs_logo](documentation/docs/images/tapirs_seq.png)
# Tapirs

Tapirs is a reproducible modular workflow for the analysis of metabarcoding data.

Tapirs is compartmentalised into several modules, all contained in the rules/ directory. Each module performs a step of the workflow.

Each rule is assigned to a Conda environment containing the appropriate packages needed to perform. Using Conda ensures version control and prevents workflow failure through package incompatability.

Tapirs was created by the EvoHull group, the University of Hull, UK



### Instructions

#### CONDA and Snakemake

First, ensure a version of conda is intalled on your machine. This can be either anaconda or miniconda.
Miniconda is available at: https://docs.conda.io/en/latest/miniconda.html

Next, ensure snakemake is installed in your base conda environment.
This can be done with the following command:

$ conda activate base
$ conda install -c bioconda -c conda-forge snakemake


#### Database setup

##### BASTA
From the base conda environment, ensuring you are in the "Tapirs/" directory, run the following commands:

$ conda env create --file envs/basta_LCA.yaml
$ conda activate basta_LCA
$ basta download gb
$ basta taxonomy
$ conda deactivate

##### BLAST
BLAST requires a local custom database for the workflow to run efficiently.
######## GS or DL - can you put a little bit in here about maing blast databases etc - Mike


#### Input wrangling

Place all library directories within the "data/" directory, ensuring they have the suffix ".fastq.gz".

Open the configuration file "config.yaml" and follow steps there.

Once this file has been configured we are ready to run the workflow.


#### Running the workflow

The workflow can be ran with the following line of code:

$ snakemake -s snakefile --use-conda --printshellcmds -n

The above line will start a dryrun. Assess the joblist of the dryrun and the fnal output files to ensure that your libraries and samples are being detected by the iteration steps.

If all looks well, run the workflow for real with the following code:

$ snakemake -s snakefile --use-conda --printshellcmds

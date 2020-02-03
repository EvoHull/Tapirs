
!!! note
    These documents assume a unix system like OSX or Linux

You will need some very basic knowledge of the command line and Snakemake. Snakemake is a relatively easy workflow manager, but we recommend that you familiarise yourself with it, perhaps carry out the [tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

# install miniconda

Conda (miniconda) is a package manager and is required to install software and their dependencies. Follow the [installation instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) for Miniconda for your operating system.

# create and activate a conda environment

Navigate to the tapirs directory.

`conda create env -f environment.yaml`

`conda activate tapirs`

The software listed in environment.yaml should now be installed. You also need to install BASTA, KronaTools, and a series of databases for blastn, BASTA, Krona, and Kraken2. Databases are large files, and everyone needs a different one, so they are not included with this install. BASTA can be problematic to install and requires its own python environment separate from the rest of the workflow.

# create a second environment, for BASTA

!!! warning
    do we need to do this or will snakemake take care of creating the environment?

The LCA analysis software BASTA requires a different version of python (2.7) and needs to be run in its own environment. Snakemake will take of this but you should first create it from the `basta_env.yaml` file.

`conda create env -f envs/basta_env.yaml`

# git clone tapirs

!!! warning
    do we need git in the environment file?


clone the repository to your local environment with `git clone github_address`

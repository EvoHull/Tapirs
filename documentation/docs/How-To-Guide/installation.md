
!!! note
    These documents assume a unix system like OSX or Linux

Although you can install and run Tapirs without too many steps you will need some very basic knowledge of the command line. A basic knowledge of Snakemake will help you to modify and configure Tapirs. Snakemake is a relatively easy workflow manager, but we recommend that you familiarise yourself with it, perhaps carry out the [tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

# git clone tapirs
If you have git installed then clone the repository to your local environment with `git clone github_address`

You could alternatively download the repository from the [Tapirs github repository](https://github.com/davelunt/Tapirs) using the green button. Then expand the zip file and navigate into the directory.

# install miniconda

Conda (miniconda) is a package manager and is required here to install software and their dependencies. Follow the [installation instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) for Miniconda for your operating system.

# create and activate a conda environment

Navigate to the tapirs directory.

`conda create env -f environment.yaml`

`conda activate tapirs`

The software listed in environment.yaml file should now be installed. You also need to install BASTA, KronaTools, and a series of databases for blastn, BASTA, Krona, and Kraken2. Databases are large files, and everyone needs a different one, so they are not included with this install. BASTA can be problematic to install and requires its own python environment separate from the rest of the workflow.

# create a second environment, for BASTA

!!! warning
    do we need to do this or will snakemake take care of creating the environment?

The LCA analysis software BASTA requires a different version of python (2.7) and needs to be run in its own environment. Snakemake will take of this but you should first create it from the `basta_env.yaml` file.

`conda create env -f envs/basta_env.yaml`

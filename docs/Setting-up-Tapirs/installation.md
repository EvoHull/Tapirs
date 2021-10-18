# INSTALLATION

!!! note
    These documents assume a unix system like OSX or Linux, on MS Windows it may be best to work in a Linux Virtual Machine

Although you can install and run Tapirs without too many steps you will need some very basic knowledge of the command line. A basic knowledge of Snakemake will help you to modify and configure Tapirs. Snakemake is a relatively easy workflow manager, but we recommend that you familiarise yourself with it, perhaps carry out the [tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

## Install miniconda

Conda (miniconda, Anaconda) is a package and environment manager and is required here to install software and their dependencies. If you choose to install Anaconda (instead of miniconda) it will also install a lot of scientific software packages. Miniconda is much more lightweight and is our recommended option.

Follow the [installation instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) for Miniconda for your operating system.

- If unsure on OSX choose "Miniconda3 MacOSX 64-bit pkg" (or similar name ending in pkg) as this gives a typical Apple install gui.

## git clone Tapirs

Apple OSX and Linux should both come with git already installed. At the command line type `git --version` and you should see the version number. If instead it reports `command not found: git` or similar then it is not installed. You can go to the [git website](https://git-scm.com/) to get installation advice or slightly easier might be to try `conda install git` at the command line.

If you have git installed then clone the repository to your local machine with `git clone https://github.com/EvoHull/Tapirs.git`

You could alternatively download the repository from the [Tapirs github repository](https://github.com/EvoHull/Tapirs) using the green "Clone or download" button. Then expand the zip file and navigate into the directory.

## Create a working environment

It is best practice to install the software you require for a specific software project (eg Tapirs) in a dedicated 'environment'. This environment will contain only the software you choose for this project and hopefully avoid software conflicts. Conda is the tool for creating and using software environments.

When snakemake runs it can be told to create separate environments to optimally run each rule (a "sequence-quality-control" environment, a "blast" environment, a "kraken2" environment). Currently however the software requirements for Tapirs are relatively straightforward and we have everything in a single "tapirs" environment. The software required in this environment is specified in the `workflow/envs/env.yaml` file (we also have an alternative tapirs.yaml for development). Clearly specifying software in this way is an important component of reproducibility.

The first time you run Tapirs it can take a while (>10 minutes) to download software and create the environment, subsequent runs will not require this step. We recommend that you create and populate this environment before you start running Tapirs.

### Install all software from the environment.yaml list

Although you can install software packages one at a time it is not very efficient. Instead we have created a list of the required software in a text file in the `workflow/envs/` directory called `env.yaml`. Conda can be told to create or update your environment by reading all the software packages this file.

Install all required software now:

`conda env create --file workflows/envs/env.yaml`

This will take a few minutes to install all the required software and their dpendencies. If in future you modify the environment file you can update it with:

`conda env update -f workflows/envs/env.yaml`

Accept the defaults (Yes) of any install questions you afre asked.

You will need to make sure that the 'tapirs' environment is active, or else the required software will not be available to be run by snakemake.

`conda activate tapirs`

If you are unsure what environments you have, and which is active, you can run:

`conda info --envs`

If you get errors when running Tapirs suggesting that some "software-name" is unknown it is most likely an issue with environments, start by checking that the "tapirs" environment is active.

The Snakemake workflow manager software was listed in the `env.yaml` file and has already been installed if you have carried out the instructions above. You could test this with a `snakemake -help` command. If you get an error such as `command not found: snakemake` its likely that the tapirs environment is not active, try: `conda activate tapirs`

## Testing the Installation

You should now have installed all the software required for your analysis. You will also require some data and reference databases. Below we give instructions of running Tapirs with test data and databases in oerder to verify the installation. Subsequently you will remove the outputs of the test and repeat with your own data and databases.

### Install the taxonomy data
The three commands below will download and expand the taxonomy data

`wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip`
`unzip new_taxdump.zip -d resources/databases/new_taxdump`
`rm new_taxdump.zip`

### Dry run Tapirs

We have provided test data
- a list of samples `/config/hull_test.tsv`
- blast and kraken2 databases in `resources/databases/`
- sequence data in `resources/libraries`

Dry run Tapirs with:

`snakemake -npr`

If you don't get any obvious errors but rather a (yellow) list of jobs to be carried out then you should be ready to go.

`snakemake --cores 4`

You can omit the number to run it on all availabale cores. This will begin the analysis of the test data. It will only take a couple of minutes to complete. Again examine the output for obvious errors (usually in red text) but everything is likely fine if it completes with something like:
```
Finished job 0.
29 of 29 steps (100%) done
```

### remove test output and prepare to run your own data
The output from any snakemake run can be cleaned up with `--delete-alloutput` but we recommend a dry run first to see what it will delete

`snakemake --delete-all-output -n`

It should list files in the `results` directory, and you can delete them all with 

`snakemake --delete-all-output`


Instructions of what files are required are provided on the Tapirs [setup](setup.md) page.

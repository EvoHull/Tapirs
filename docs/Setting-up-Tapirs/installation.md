
!!! note
    These documents assume a unix system like OSX or Linux, but should be applicable to all systems

Although you can install and run Tapirs without too many steps you will need some very basic knowledge of the command line. A basic knowledge of Snakemake will help you to modify and configure Tapirs. Snakemake is a relatively easy workflow manager, but we recommend that you familiarise yourself with it, perhaps carry out the [tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

# Install miniconda

Conda (miniconda, Anaconda) is a package and environment manager and is required here to install software and their dependencies. If you choose to install Anaconda (instead of miniconda) it will also install a lot of scientific software packages. Miniconda is much more lightweight and is our recommended option.

Follow the [installation instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) for Miniconda for your operating system.

- If unsure on OSX choose "Miniconda3 MacOSX 64-bit pkg" (or similar name ending in pkg) as this gives a typical Apple install gui.

# git clone tapirs
Apple OSX and Linux should both come with git already installed. At the command line type `git --version` and you should see the version number. If instead it reports `command not found: git` or similar then it is not installed. You can go to the [git website](https://git-scm.com/) to get installation advice or slightly easier might be to try `conda install git` at the command line.

If you have git installed then clone the repository to your local machine with `git clone https://github.com/EvoHull/Tapirs.git`

You could alternatively download the repository from the [Tapirs github repository](https://github.com/EvoHull/Tapirs) using the green "Clone or download" button. Then expand the zip file and navigate into the directory.

# Create a working environment called 'tapirs'
It is best practice to install the software you require for a specific software project (eg Tapirs) in a dedicated 'environment'. This environment will contain only the software you choose for this project and hopefully avoid software conflicts. Conda is the tool for creating and using software environments.
```
conda create --name tapirs
conda activate tapirs
```

When snakemake runs it can be told to create separate environments to optimally run each rule (a "sequence-quality-control" environment, a `blast` environment, a `kraken2` environment). Currently however the software requirements for Tapirs are relatively straightforward and we have everything in a single 'tapirs' environment. The software required in this environment is specified in the `envs/environment.yaml` file. 

The first time you run Tapirs it can take a while (>10 minutes) to download software and create the environment, subsequent runs will not require this step. We recommend that you create and populate this environment before you start running Tapirs.

## Install all software from the environment.yaml list
Although you can install software packages one at a time it is not very efficient. Instead we have created a list of the required software in a text file in the `envs/` directory. Conda can be told to update your environment by reading all the software packages this file.

If you have not yet told conda to create a 'tapirs' environment then you can do so now, and install all required software:

`conda env create --file envs/environment.yaml`

If you have, above, created a 'tapirs' environment you can update it to include all the required software with:

`conda env update -f environment.yaml`

Accept the defaults (Yes) of any install questions you afre asked.

You will need to make sure that the 'tapirs' environment is active, or else the required software will not be available to be run by snakemake.

`conda activate tapirs`

If you are unsure what environments you have, and which is active, you can run:

`conda info --envs`

If you get errors when running Tapirs suggesting that some software-name is unknown it is most likely an issue with environments, start by checking that the 'tapirs' environment is active.

The Snakemake workflow manager software was listed in the `environment.yaml` file and has already been installed if you have carried out the instructions above. You could test this with a `snakemake -help` command. If you get an error such as `command not found: snakemake` its likely that the tapirs environment is not active, try: `conda activate tapirs`

## Optional: install Snakemake if you wish snakemake to handle your envs

If you have decided to let snakemake itself handle all software environments as it runs, and have not installed the 'tapirs' environment above from the `environment.yaml` file, you will need to install snakemake now.

```
conda install -c conda-forge snakemake
```
test your install with `snakemake -help`

This approach will work well, and in some instances will be preferable, but:
1. You must remember to specify conda when running snakemake, e.g.
`snakemake --use-conda --cores 8`
2. You will need to be patient the *first time* you run a command to use snakemake as it must create the environments and downlaod their software before beginning the analysis workflow.

# Databases and data
You should now have installed all the software required for your analysis. You will also require some data and databases however.

In order to search a reference database with your query sequences you will need to provide both, and tell Tapirs where they are located. Databases are large files, and everyone needs a different one, so they are not included with Tapirs.

Instructions of what files are required are provided on the Tapirs [setup](setup.md) page.

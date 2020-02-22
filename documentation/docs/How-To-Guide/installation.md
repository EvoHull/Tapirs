
!!! note
    These documents assume a unix system like OSX or Linux

Although you can install and run Tapirs without too many steps you will need some very basic knowledge of the command line. A basic knowledge of Snakemake will help you to modify and configure Tapirs. Snakemake is a relatively easy workflow manager, but we recommend that you familiarise yourself with it, perhaps carry out the [tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

# install miniconda

Conda (miniconda) is a package manager and is required here to install software and their dependencies. Follow the [installation instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) for Miniconda for your operating system.

# git clone tapirs
Make sure that you have git installed. At the command line type `git --version` and you should see the version number. If instead it reports `command not found: git` or similar then it is not installed. You can go to the [git site](https://git-scm.com/) to get installation advice or slightly easier might be to try `conda install git` at the command line.

If you have git installed then clone the repository to your local environment with `git clone https://github.com/davelunt/Tapirs.git`

You could alternatively download the repository from the [Tapirs github repository](https://github.com/davelunt/Tapirs) using the green button. Then expand the zip file and navigate into the directory.

# create and activate a conda environment
Make sure you are in the tapirs directory, then give these 2 commands:

```
conda create env -f envs/tapirs.yaml
conda activate tapirs
```

The software listed in the `tapirs.yaml` environment file should now be installed.

# KronaTools install
Unfortunately you will need to complete the installation of KronaTools semi-manually with an install script. [Advice on installing Krona Tools](https://github.com/marbl/Krona/wiki/Installing)] is give on their Github wiki.

Instructions are also given on the Tapirs [setup page](setup.md). Make sure that you have issued the two commands:

```
mkdir data/databases/krona/
ktUpdateTaxonomy.sh data/databases/krona/
```

# Databases and data
In order to search a reference database with your query sequences you will need to provide both. Databases are large files, and everyone needs a different one, so they are not included with Tapirs.

Instructions of what files are required are provided on the Tapirs [setup](setup.md) page.

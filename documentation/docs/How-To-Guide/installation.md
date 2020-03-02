
!!! note
    These documents assume a unix system like OSX or Linux

Although you can install and run Tapirs without too many steps you will need some very basic knowledge of the command line. A basic knowledge of Snakemake will help you to modify and configure Tapirs. Snakemake is a relatively easy workflow manager, but we recommend that you familiarise yourself with it, perhaps carry out the [tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

# Install miniconda

Conda (miniconda) is a package manager and is required here to install software and their dependencies. Follow the [installation instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) for Miniconda for your operating system.

# git clone tapirs
Apple OSX and Linux should both come with git already installed. At the command line type `git --version` and you should see the version number. If instead it reports `command not found: git` or similar then it is not installed. You can go to the [git site](https://git-scm.com/) to get installation advice or slightly easier might be to try `conda install git` at the command line while within your conda base environment.

If you have git installed then clone the repository to your local machine with `git clone https://github.com/davelunt/Tapirs.git`

You could alternatively download the repository from the [Tapirs github repository](https://github.com/davelunt/Tapirs) using the green "Clone or download" button. Then expand the zip file and navigate into the directory.

# Create and activate a conda environment
Make sure you are in the tapirs directory, then give these 2 commands:

```
conda create env -f envs/tapirs.yaml
conda activate tapirs
```

The first command has created a new environment (called "tapirs") containing all the software specified in the tapirs.yaml file. The second command has activated that new environment, making all the software available to you.

# KronaTools install
Krona is used to produce active html reports on the taxonomy. You will need to complete the installation of KronaTools semi-manually with an install script. [Advice on installing Krona Tools](https://github.com/marbl/Krona/wiki/Installing) is given on their Github wiki.

Instructions are also given on the Tapirs [setup page](setup.md). Make sure that you have issued the two commands:

```
mkdir data/databases/krona/
ktUpdateTaxonomy.sh data/databases/krona/
```

# Databases and data
You should now have installed all the software required for your analysis. You will also require some data and databases however.

In order to search a reference database with your query sequences you will need to provide both, and tell Tapirs where they are located. Databases are large files, and everyone needs a different one, so they are not included with Tapirs.

Instructions of what files are required are provided on the Tapirs [setup](setup.md) page.

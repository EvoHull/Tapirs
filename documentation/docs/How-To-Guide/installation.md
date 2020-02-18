
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

Navigate to the tapirs directory.

`conda create env -f environment.yaml`

`conda activate tapirs`

The software listed in environment.yaml file should now be installed.

# Additional databases and installs
Unfortunately you will need to complete the installation of KronaTools semi-manually with an install script. [Advice on installing Krona Tools(https://github.com/marbl/Krona/wiki/Installing)] is give on their Github wiki.

In order to search a blast database, or a Kraken database, you will need to provide those databases yourself. Databases are large files, and everyone needs a different one, so they are not included with this install or Tapirs.

Your databases should be placed into the `data/databases/blast` `data/databases/kraken` directories.

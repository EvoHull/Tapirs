
!!! note
    These documents assume a unix system like OSX or Linux, but should be applicable to all systems

Although you can install and run Tapirs without too many steps you will need some very basic knowledge of the command line. A basic knowledge of Snakemake will help you to modify and configure Tapirs. Snakemake is a relatively easy workflow manager, but we recommend that you familiarise yourself with it, perhaps carry out the [tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

# Install miniconda

Conda (miniconda) is a package and environment manager and is required here to install software and their dependencies. Follow the [installation instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) for Miniconda for your operating system.

# git clone tapirs
Apple OSX and Linux should both come with git already installed. At the command line type `git --version` and you should see the version number. If instead it reports `command not found: git` or similar then it is not installed. You can go to the [git site](https://git-scm.com/) to get installation advice or slightly easier might be to try `conda install git` at the command line while within your conda base environment.

If you have git installed then clone the repository to your local machine with `git clone https://github.com/EvoHull/Tapirs.git`

You could alternatively download the repository from the [Tapirs github repository](https://github.com/EvoHull/Tapirs) using the green "Clone or download" button. Then expand the zip file and navigate into the directory.

Another way would to be to download the GitHub Desktop application and proceed from there.

# Install Snakemake in your base conda environment
Once conda is installed you will see "(base)" at the start of your command prompt. If this is not the case run `conda activate base`. With (base) activated, install Snakemake:

```
conda install -c bioconda -c conda-forge snakemake
```
Snakemake is the only dependency that requires manually installing; it will install and manage all other packages itself using conda.

# Databases and data
You should now have installed all the software required for your analysis. You will also require some data and databases however.

In order to search a reference database with your query sequences you will need to provide both, and tell Tapirs where they are located. Databases are large files, and everyone needs a different one, so they are not included with Tapirs.

Instructions of what files are required are provided on the Tapirs [setup](setup.md) page.

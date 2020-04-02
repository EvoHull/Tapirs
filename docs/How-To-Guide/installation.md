
!!! note
    These documents assume a unix system like OSX or Linux, but should be applicable to all systems

Although you can install and run Tapirs without too many steps you will need some very basic knowledge of the command line. A basic knowledge of Snakemake will help you to modify and configure Tapirs. Snakemake is a relatively easy workflow manager, but we recommend that you familiarise yourself with it, perhaps carry out the [tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

# Install miniconda

Conda (miniconda, Anaconda) is a package and environment manager and is required here to install software and their dependencies. If you choose to install Anaconda (instead of miniconda) it will also install a lot of scientific software packages. Miniconda is much more lightweight and is our recommended option.

Follow the [installation instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) for Miniconda for your operating system.

- If unsure on OSX choose "Miniconda3 MacOSX 64-bit pkg" (or similar name ending in pkg) as this gives a typical Apple install gui.

# git clone tapirs
Apple OSX and Linux should both come with git already installed. At the command line type `git --version` and you should see the version number. If instead it reports `command not found: git` or similar then it is not installed. You can go to the [git site](https://git-scm.com/) to get installation advice or slightly easier might be to try `conda install git` at the command line.

If you have git installed then clone the repository to your local machine with `git clone https://github.com/EvoHull/Tapirs.git`

You could alternatively download the repository from the [Tapirs github repository](https://github.com/EvoHull/Tapirs) using the green "Clone or download" button. Then expand the zip file and navigate into the directory.

# Create a working environment called 'tapirs'
It is best practice to install the software you require for a specific software project (eg Tapirs) in a dedicated 'environment'. This environment will contain only the software you choose for this project and hopefully avoid software conflicts. Conda is the tool for creating and using software environments.
```
conda create --name tapirs
conda activate tapirs
```
# Install Snakemake
You are now in an environment (called tapirs) but need to next install the software that will be used by Tapirs workflows. Snakemake can do this for you once it is installed. Install Snakemake in your conda environment:

```
conda install -c conda-forge snakemake
```
test your install with `snakemake --help`

# Install all software from the environment.yaml list
Although you can install software packages one at a time, as you just did for snakemake, it is not very efficient. Instead we have created a list of the required software in a text file in the `envs/` directory. Conda can be told to update your environment by reading all the software packages this file.
```
conda env update -f envs/environment.yaml
```
Accept defaults by answering Yes to any questions. This installation might take a few mins to complete (unless conda or the internet is very slow that day).

# Databases and data
You should now have installed all the software required for your analysis. You will also require some data and databases however.

In order to search a reference database with your query sequences you will need to provide both, and tell Tapirs where they are located. Databases are large files, and everyone needs a different one, so they are not included with Tapirs.

Instructions of what files are required are provided on the Tapirs [setup](setup.md) page.

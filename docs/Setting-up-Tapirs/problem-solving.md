This file outlines problems that users have faced and solved. Please contribute more by editing this document or raising a github issue. We cannot cover everything here of course and assume a basic understanding.

## conda build environment fails
This is rare but it is possible that one of the softwares in the environment.yaml file will fail to install. This is an issue with conda not Tapirs, try installing the program in question at the command line with `conda -c conda-forge install PROGRAM`. It might help to try a different repository such as `-c bioconda` instead of conda-forge. It might be worth googling "conda PROGRAM" just to see how Anaconda Cloud page suggests installing.

## PROGRAM not found
Please check that you are in the correct conda environment. `conda activate tapirs` should make sure that you are in the tapirs environment.

## error involving cores
Some hosts require you to specify the number of cores on which to run your job. You can just modify your snakemake command to include `--use-cores 8` although you can vary this number.

## I really want to speed up the first run building the environments
We recommend using [mamba](https://github.com/mamba-org/mamba) not conda to build environments. Snakemake fully supports mamba. It is a LOT faster.

### pre-build the environments yourself
Currently Tapirs uses a single environment file but this might change soon. If the environment is specified in a **single file** (eg `envs/environment.yaml`) you can build the environment before first run with the command `mamba env create -f env/environment.yaml`

If the environments are in multiple files (eg `envs/blast.yaml`, `envs/qc.yaml`) then these can all be installed at once from the command line `for f in envs/*.yaml; do mamba env create -f $f; done`

### pre-build the environments in snakemake
Snakemake can pre-build the environments for you, and can be told to use mamba as it's default. **These commands should work, but require full testing**

```
snakemake -conda-frontend mamba
snakemake –use-conda -conda-create-envs-only True
```
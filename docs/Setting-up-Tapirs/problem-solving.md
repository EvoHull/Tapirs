This file outlines problems that users have faced and solved. Please contribute more by editing this document or raising a github issue. We cannot cover everything here of course and assume a basic understanding.

## conda build environment fails
This is rare but it is possible that one of the softwares in the environment.yaml file will fail to install. This is an issue with conda not Tapirs, try installing the program in question at the command line with `conda -c conda-forge install PROGRAM`. It might help to try a different repository such as `-c bioconda` instead of conda-forge. It might be worth googling "conda PROGRAM" just to see how Anaconda Cloud page suggests installing.

## PROGRAM not found
Please check that you are in the correct conda environment. `conda activate tapirs` should make sure that you are in the tapirs environment.

## error involving cores
Some hosts require you to specify the number of cores on which to run your job. You can just modify your snakemake command to include `--use-cores 8` although you can vary this number.
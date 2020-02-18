
!!! warning "install first then setup"
    These instructions are to set up Tapirs after [installation](installation.md) has been carried out.

# CONFIG FILES
The `config.yaml` file contains settings for the workflow's operation. In an ideal experiment you would never have to alter any of the snakemake workflow code, only set up a configuration text file. Most of the settings in `config.yaml` have reasonable default values, but some may require your input. Open `config.yaml` in a text editor.

1. set a name for your experiment (default=expt_name)
2. set the directory containing your **demultiplexed** fastq.gz data (default=data/demultiplexed)
3. check the location of your databases is correct and amend if they are elsewhere. Defaults are:
    - `data/databases/blast`
    - `data/databases/kraken`
    - `data/databases/taxonomy`
4. set fastp parameters for quality control and read merging
5. set blast parameters
6. set MLCA parameters

# DATA

The `/data` directory should contain 2 subdirectories `/demultiplexed` `/databases`

## demultiplexed
The start point of this workflows is demultiplexed fastq.gz sequence files in the data/demultiplexed directory. There are many ways to produce these files, and your sequencing machine or sequencing centre will most likely return this as the default data.

It is essential for reproducibility that this starting data is kept together with the experiment. If size prohibits its inclusion then archive it publicly and include the doi into the experimental record.

!!! tip "tip: exclude data directory from version control"
    Remember that if you are wisely keeping your analysis under git version control, then you can exclude very large data files from syncing to the web by including the data directory in your .gitignore file. It is essential however that you make other arrangements for both backing up and publishing the data.

## Databases
### Kraken2
In order to run a Kraken analysis you will need to create a Kraken database for the program to search your query sequences against. These databases can be large and require a lot of RAM and time to build them. They only need to be made once however and then all subsequent analyses can be performed against the same database. If you are building a custom database from a few thousand sequences then database construction will likely be quick.

See the [Kraken Tutorial](../Tutorials/kraken_tutorial.md) and the [Kraken2 manual](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual) for information on building Kraken databases.

!!! note "create a Kraken database"
    If you have a single fasta format file (allseqs.fas) containing all the sequences to be included in the database, then you could create a kraken database with the command:

    `database creation example here`

It is essential for reproducibility that you publicly archive your database at the end of the experiment. Zenodo.org is a suitable location.

### blast
You will need to build a blastn database from a collection of fasta files. Information on this can be found at the NCBI site [Blast help pages](https://www.ncbi.nlm.nih.gov/books/NBK279680/).

!!! note "create a blast database"
    If you have a single fasta format file (allseqs.fas) containing all the sequences to be included in the database, then you could create a blast database with the command:

    `database creation example here`
See the [Blast Tutorial](../Tutorials/blast_tutorial.md) for more detailed help

It is essential for reproducibility that you publicly archive your database at the end of the experiment. [Zenodo.org](zenodo.org) is a suitable location.

### taxonomy database
Several programs require the NCBI taxonomy database in order to carry out taxonomic assignment (BASTA, Kraken).

# DRY RUN SNAKEFILE
Make sure you are in the directory containing the snakefile then type `snakefile -npr` to dry-run the workflow.

If all has gone well Snakemake will report the jobs it needs to perform without any complaint. If not (and it is the usual situation) you will need to diagnose and fix any minor issues.

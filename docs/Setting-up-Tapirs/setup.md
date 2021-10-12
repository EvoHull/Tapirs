# SETUP

!!! warning "install first then setup"
    These instructions are to set up Tapirs for your experiment after [installation](installation.md) has been carried out.

## OVERVIEW

It is vital to carry out all the set up instructions. This can be slow, and require some thought, but it only has to be done once.

1. edit the `config.yaml` file for this experiment, it controls the workflow
2. make sure your data is in the right location and format
3. create the databases against which your query sequences will be compared.

## CONFIG FILES

The `config/config.yaml` file contains settings for the workflow's operation. In an ideal experiment you would never have to alter any of the snakemake workflow code, only set up a configuration text file. Most of the settings in `config.yaml` have reasonable default values, but some may require your input. Open `config.yaml` in a text editor.

1. set a name for your experiment
2. check that the location of your raw data is correctly specified (`resources/libraries/LIBRARY/SAMPLE.fastq.gz`)
3. check the locations of your reference sequence databases are correctly specified, for example
    - `resources/databases/blast`
    - `resources/databases/kraken2`
    - `resources/databases/new_taxdump`

In addition you may wish to fine-tune the analysis parameters of the programs. These parameters have reasonable defaults and changing them is not compulsory (unlike the options above, where you _must_ identify your data).

- set fastp parameters for quality control and read merging
- set blast parameters
- set MLCA parameters

## DATA

The `resources/` directory must contain 2 subdirectories: a) `libraries/` b) `databases/`

## demultiplexed sequence data

The start point of this workflow is demultiplexed `fastq.gz` sequence files in the `resources/libraries` directory. Each directory in `libraries` is a sequencing library (eg `resources/libraries/river1`) and each library contains your sequences in 1 or more fastq files (eg `sample1.fastq.gz`)

It is essential for reproducibility that this starting data is kept together with the experiment. If size prohibits its inclusion then archive it publicly and include the doi into the experimental record.

!!! tip "tip: exclude data directory from version control"
    Remember that if you are wisely keeping your analysis under git version control, then you can exclude very large data files from syncing to the web by including the data directory in your .gitignore file. It is essential however that you make other arrangements for both backing up and publishing the data.

To facilitate iteration, the input files must follow the structure of `yourlibrary/yoursample.R*.fastq.gz`, the full path being `Tapirs/resources/libraries/yourlibrary/yoursample.R*.fastq.gz`. The * asterisk stands for the read number, ie R1 or R2.

### libraries and samples lists

You must make sure that your config file points to your list of samples and libraries (eg samples: `config/hull1.tsv`)

Rather than implement a complex set of wildcards the recommended approach for snakemake is to create textual lists (.tsv files) of the sequencing metadata [see "samples, libraries, and units .tsv" Note below]. These lists are generated for you by a script, but this approach allows quality control and human intervention in determining what samples are analysed in each run. Importantly having an explicit textual record of what samples were analysed, and their metadata, is good for reproducibility.

!!! Note "Note: samples, libraries, and units .tsv"
    The naming of the files differs between biological disciplines. In our work a library often represents a physical location ("Lake Windermere") and a sample represents a physical unit taken for DNA extraction ("water-sample-12"). In other disciplines the meaning of "sample" may differ slightly.

Tapirs `workflow/scripts` contains a python script `get_dirs_files.py` to create a text file `samples.tsv`. This contains a list of all libraries and sample names. You can hand edit and should sanity check this list of input files, then specify it in the `config.yaml`. Run the python script at the command line with `python workflow/scripts/get_dirs_files.py` to generate the tsv file

## DATABASES

### Kraken2

In order to run a Kraken analysis you will need to create a Kraken database for the program to search your query sequences against. These databases can be large and require a lot of RAM and time to build them. They only need to be made once however and then all subsequent analyses can be performed against the same database. If you are building a custom database from a few thousand sequences then database construction will likely be quick.

Creating a local Kraken2 database for Tapirs requires a fasta input file containing all reference sequences.

Database creation has 3 steps:

1. Download taxonomy:
`kraken2-build --download-taxonomy --db fish_kraken --threads 6`
2. Add sequences to library:
`kraken2-build --add-to-library databases/12s_verts.fasta --no-masking --db fish_kraken --threads 6`
the `--no-masking` flag improved accuracy for short reads
3. Build database:
`kraken2-build --build --minimizer-spaces 1 --kmer-len 31 --db fish_kraken --threads 6`

it is important to have some minimizer-spaces in lmers, having 1 makes for more accurate taxonomic assignment. Having kmer same length as lmer for shorter reads makes for more accurate taxonomic assignment

`fish_kraken` and `12s_verts.fasta` are our experiment-specific names for our data, you will choose different names for whatever data you are working with.

See the [Kraken Tutorial](../Tutorials/kraken_tutorial.md) and the [Kraken2 manual](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual) for information on building Kraken databases.

It is essential for reproducibility that you publicly archive your database at the end of the experiment. [Zenodo.org](https://zenodo.org) is a suitable location. If your database is too large to database easily, consider making it from scripts in a reproducible manner instead and then archiving those scripts. Don't forget to include the international sequence database **version** from which your scripts harvested the sequences.

### blast

You will need to build a blastn database from a collection of fasta files. Information on this can be found at the NCBI site [Blast help pages](https://www.ncbi.nlm.nih.gov/books/NBK279680/).

You will require:

- Fasta input file containing all reference sequences
- Accession to taxid map. See [NCBI blast instructions](https://www.ncbi.nlm.nih.gov/books/NBK279688/) for more details.

!!! note "create a blast database"
    If you have a single fasta format file (allseqs.fas) containing all the sequences to be included in the database, then you could create a blast database with the command:

    ```
    makeblastdb -in blast_db/allseqs.fasta -input_type fasta -dbtype nucl \
    -parse_seqids -taxid_map blast_db/tax_map.txt -out blast_db/allseqs
    ```

See the [Blast Tutorial](../Tutorials/blast_tutorial.md) for more detailed help

It is essential for reproducibility that you publicly archive your database at the end of the experiment. [Zenodo.org](https://zenodo.org) is a suitable location.

### taxonomy database

Several programs may require the NCBI taxonomy database in order to carry out taxonomic assignment. This is typically contained in a `new_taxdump` directory that is fetched directly from the NCBI.

`mkdir resources/databases/new_taxdump`

`cd resources/databases/new_taxdump`

`wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip`

`unzip new_taxdump`

Alternatively there is a snakemake rule specifically to set this up. You can run this with the command:
`snakemake -s workflow/rules/tax_db.smk`

If on OSX `wget` is not installed, you can install with `conda install -c conda-forge wget` or any ftp client will fetch this new_taxdump data for you.

## DRY RUN TAPIRS

Make sure you are in the directory containing the snakefile then type `snakemake -npr` to dry-run the workflow.

If all has gone well Snakemake will report the jobs it needs to perform (in yellow) without any complaint. If not (as is common in most experiments) you will need to diagnose and fix any minor issues. Some errors are only detected in the real run, not the dry run, and they often concern the format of data files, as these have not been checked by a dry run.

## RUN TAPIRS

Run Tapirs with `snakemake --cores 4` or `snakemake --cores 4 --printshellcmds` commands

Tapirs should now run, processing the data, assigning taxonomy using blast and/or kraken2, and writing reports.

When it finishes you should also ask it to write a report with the command
`snakemake --report reports/snakemake_report.html`

## EXCLUDE ANALYSES

If you wish to run Tapirs without invoking one of analysis programs (eg Kraken2 or blast) then you can specify this in the config file. For example replace the line `analysis_method: "both"` with `analysis_method: "blast"` to restrict the analysis to just blast.

Remember that snakemake does not rerun jobs already completed. So if you run just blast, then run just kraken2, it will not attempt to repeat the qc stages common to both unless you change the input data.

## REMOVING FILES FROM PREVIOUS RUNS

Snakemake can clean up files it as previously created. This is useful if you have reports and intermediate results from previous runs that you wish to remove before a new run. The Snakemake docs have a [FAQ on cleaning files](https://snakemake.readthedocs.io/en/stable/project_info/faq.html#how-do-i-remove-all-files-created-by-snakemake-i-e-like-make-clean), in short though try `snakemake some_target --delete-all-output --dry-run` The`--dry-run` flag checks what will be removed before you do it, when it looks fine rerun without `--dry-run`.

We highly recommend performing a `--dry-run` as `--delete-all-output` is as dangerous to your results as it sounds.

## REFERENCES

Altschul, S. F. et al. (1990) ‘Basic local alignment search tool’, Journal of molecular biology, 215(3), pp. 403–410. [doi: 10.1016/S0022-2836(05)80360-2](https://doi.org/10.1016/S0022-2836(05)80360-2)

Chen, S. et al. (2018) ‘fastp: an ultra-fast all-in-one FASTQ preprocessor’, Bioinformatics . Oxford University Press, 34(17), pp. i884–i890. [doi: 10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)

Daniel McDonald, Jose C. Clemente, Justin Kuczynski, Jai Ram Rideout, Jesse Stombaugh, Doug Wendel, Andreas Wilke, Susan Huse, John Hufnagle, Folker Meyer, Rob Knight, and J. Gregory Caporaso. The Biological Observation Matrix (BIOM) format or: how I learned to stop worrying and love the ome-ome. GigaScience 2012, 1:7. [doi:10.1186/2047-217X-1-7](https://doi.org/10.1186/2047-217X-1-7)

Ondov, B. D., Bergman, N. H. and Phillippy, A. M. (2011) ‘Interactive metagenomic visualization in a Web browser’, BMC bioinformatics, 12, p. 385. [doi: 10.1186/1471-2105-12-385](https://doi.org/10.1186/1471-2105-12-385)

Rognes, T. et al. (2016) ‘VSEARCH: a versatile open source tool for metagenomics’, PeerJ, 4, p. e2584. [doi: 10.7717/peerj.2584]

Wood, D. E., Lu, J. and Langmead, B. (2019) ‘Improved metagenomic analysis with Kraken 2’, Genome biology, 20(1), p. 257. [doi: 10.1186/s13059-019-1891-0](https://doi.org/10.1186/s13059-019-1891-0)

# RUNNING TAPIRS

!!! warning "install and setup before running"
    These instructions are to run Tapirs after [installation](../Setting-up-Tapirs/installation.md) and [set up](../Setting-up-Tapirs/setup.md) have been carried out.

## DRY RUN TAPIRS

Make sure you are in the top-level directory containing the snakefile then type `snakemake -npr`  or `snakemake -s snakefile --printshellcmds -n -k` to dry-run the workflow.

If all has gone well Snakemake will report the jobs it needs to perform without any complaint. If not (as is common in most experiments) you will need to diagnose and fix any minor issues. Reading the [problem-solving](../Setting-up-Tapirs/problem-solving.md) documentation might help. Some errors are only detected in the real run, not the dry run, and they often concern the format of data files, as these have not been checked by a dry run.

## RUN TAPIRS

Run Tapirs with either the `snakemake --printshellcmds --cores 4` command. Cores can be set to any number depending on your machine, and even basic machines should have 4.

Tapirs should now run, processing the data, assigning taxonomy using blast and kraken2, and writing reports.

When it finishes you should also ask it to write a report with the command

`snakemake --report reports/snakemake_report.html`

## EXCLUDE ANALYSES

If you wish to run Tapirs without invoking one of analysis programs (eg Kraken2 or blast) then you can specify this in the config file.

## REMOVING FILES FROM PREVIOUS RUNS


Snakemake can clean up files it as previously created. This is useful if you have reports and intermediate results from previous runs that you wish to remove before a new run. The Snakemake docs have a [FAQ on cleaning files](https://snakemake.readthedocs.io/en/stable/project_info/faq.html#how-do-i-remove-all-files-created-by-snakemake-i-e-like-make-clean), in short though try `snakemake some_target --delete-all-output --dry-run` The`--dry-run` flag checks what will be removed before you do it, when it looks fine rerun without `--dry-run`.

We highly recommend performing a `--dry-run` as `--delete-all-output` is as dangerous to your results as it sounds.

## **REFERENCES**

Altschul, S. F. et al. (1990) ‘Basic local alignment search tool’, Journal of molecular biology, 215(3), pp. 403–410. [doi: 10.1016/S0022-2836(05)80360-2](https://doi.org/10.1016/S0022-2836(05)80360-2)

Chen, S. et al. (2018) ‘fastp: an ultra-fast all-in-one FASTQ preprocessor’, Bioinformatics . Oxford University Press, 34(17), pp. i884–i890. [doi: 10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)

Daniel McDonald, Jose C. Clemente, Justin Kuczynski, Jai Ram Rideout, Jesse Stombaugh, Doug Wendel, Andreas Wilke, Susan Huse, John Hufnagle, Folker Meyer, Rob Knight, and J. Gregory Caporaso. The Biological Observation Matrix (BIOM) format or: how I learned to stop worrying and love the ome-ome. GigaScience 2012, 1:7. [doi:10.1186/2047-217X-1-7](https://doi.org/10.1186/2047-217X-1-7)

Ondov, B. D., Bergman, N. H. and Phillippy, A. M. (2011) ‘Interactive metagenomic visualization in a Web browser’, BMC bioinformatics, 12, p. 385. [doi: 10.1186/1471-2105-12-385](https://doi.org/10.1186/1471-2105-12-385)

Rognes, T. et al. (2016) ‘VSEARCH: a versatile open source tool for metagenomics’, PeerJ, 4, p. e2584. [doi: 10.7717/peerj.2584]

Wood, D. E., Lu, J. and Langmead, B. (2019) ‘Improved metagenomic analysis with Kraken 2’, Genome biology, 20(1), p. 257. [doi: 10.1186/s13059-019-1891-0](https://doi.org/10.1186/s13059-019-1891-0)

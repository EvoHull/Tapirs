# conda and snakemake cheatsheet
Some commands to help when you're struggling to remember the syntax:

## Conda create & activate

```
conda env create --name snake-env --file environment.yaml
```
if .yaml has name option specified
```
conda env create -f environment.yaml
```
```
conda activate snake-env
conda env update -f environment.yaml
conda info --envs
conda env export > environment_out.yml
conda install -c bioconda softwarename
```

## Snakemake
```
snakemake --dag | dot -Tsvg > dag.svg
snakemake --dag results/final_fasta/*.fastq | dot -Tsvg > dag.svg
snakemake --rulegraph | dot -Tpng > rule-graph.png

snakemake -s testsnakefile.smk
snakemake --forceall
snakemake --delete-all-output -n
```
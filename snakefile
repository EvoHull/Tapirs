#------------------------------------------------------
# Tapirs
# ---------
# A reproducible metabarcoding workflow using snakemake
#------------------------------------------------------

#configfile: "config.yaml"

# Flag "$ snakemake" with "--report" to use
report: "reports/tapirs.rst"   ### Check to make sure this works and that the output is something sensible - Mike

## Need to implement something to allow intake of both fastq and fastq.gz -Mike
## This count be a rule wherein if the input is eg. .gz, it unzips,
## whereas if not it automatically continues to the next rule that takes unzipped - mike



library = "testlib"

sample,= glob_wildcards("data/01_demultiplexed/testlib/{sample}.R1.fastq.gz")

R=["R1", "R2"]

my_experiment = "testing_tapirs"

conda_envs = ["tapirs.yaml", "basta_LCA.yaml"]

#-------------------------------------------------------------------------------
# Target rules
#-------------------------------------------------------------------------------
rule all:
    input:
# results ----------------------------------------------------------------------
        expand("results/02_trimmed/{library}/{sample}.{R}.fastq.gz", library=library, sample=sample, R=R),
        expand("results/02_trimmed/{library}/{sample}.unpaired.{R}.fastq.gz", library=library, sample=sample, R=R),
        expand("results/02_trimmed/{library}/{sample}.merged.fastq.gz", library=library, sample=sample),
        expand("results/03_denoised/{library}/{sample}.fasta", library=library, sample=sample, R=R),
        expand("results/blast/{library}/{sample}_blast.out", library=library, sample=sample),
        #expand("results/blast/{library}/{sample}_blast.taxed.out", library=library, sample=sample),
        expand("results/mlca/{library}/{sample}_lca.tsv", library=library, sample=sample),
# reports ----------------------------------------------------------------------
        expand("reports/fastp/{library}/{sample}.json", library=library, sample=sample),
        expand("reports/fastp/{library}/{sample}.html", library=library, sample=sample),
        expand("reports/vsearch/{library}/{sample}.denoise.biom", library=library, sample=sample),
        expand("reports/vsearch/{library}/{sample}_fq_eestats", library=library, sample=sample),
        expand("reports/vsearch/{library}/{sample}_fq_readstats", library=library, sample=sample),
        expand("reports/archived_envs/{conda_envs}", conda_envs=conda_envs),
        expand("results/kraken/{my_experiment}.tsv", my_experiment=my_experiment),
        expand("reports/krona/kraken/{library}/{sample}.html", library=library, sample=sample),
        expand("reports/krona/mlca/{library}/{sample}.html", library=library, sample=sample)

#-----------------------------------------------------
# Rule files
#-----------------------------------------------------
# include: "rules/reports.smk",
# include: "rules/kraken.smk",
# include: "rules/blast.smk",
# include: "rules/qc.smk"

#-----------------------------------------------------
# fastp, control for sequence quality and pair reads
#-----------------------------------------------------
# maybe this should be 2 rules, trim and then merge?
rule fastp_trim_and_merge:
    # message: "Beginning fastp quality control of raw data"
    conda:
        "envs/tapirs.yaml"
    input:
        read1 = "data/01_demultiplexed/{library}/{sample}.R1.fastq.gz",
        read2 = "data/01_demultiplexed/{library}/{sample}.R2.fastq.gz"
    output:
        out1 = "results/02_trimmed/{library}/{sample}.R1.fastq.gz",
        out2 = "results/02_trimmed/{library}/{sample}.R2.fastq.gz",
        out_unpaired1 = "results/02_trimmed/{library}/{sample}.unpaired.R1.fastq.gz",
        out_unpaired2 = "results/02_trimmed/{library}/{sample}.unpaired.R2.fastq.gz",
        out_failed = "results/02_trimmed/{library}/{sample}.failed.fastq.gz",
        merged = "results/02_trimmed/{library}/{sample}.merged.fastq.gz",
        json = "reports/fastp/{library}/{sample}.json",
        html = "reports/fastp/{library}/{sample}.html"
    shell:
        "fastp \
        -i {input.read1} \
        -I {input.read2} \
        -o {output.out1} \
        -O {output.out2} \
        --unpaired1 {output.out_unpaired1} \
        --unpaired2 {output.out_unpaired2} \
        --failed_out {output.out_failed} \
        -j {output.json} \
        -h {output.html} \
        --qualified_quality_phred 30 \
        --length_required 90 \
        --cut_tail \
        --trim_front1 20 \
        --trim_front2 20 \
        --max_len1 106 \
        --max_len2 106 \
        --merge \
        --merged_out {output.merged} \
        --overlap_len_require 90 \
        --correction \
        "

#-----------------------------------------------------
# convert files from fastq to fasta
#-----------------------------------------------------
rule fastq_to_fasta:
    conda:
        "envs/tapirs.yaml"
    input:
        "results/02_trimmed/{library}/{sample}.merged.fastq.gz"
    output:
        "results/02_trimmed/{library}/{sample}.merged.fasta",
    shell:
        "vsearch \
        --fastq_filter {input} \
        --fastaout {output} \
        "

#-----------------------------------------------------
# vsearch, fastq report
#-----------------------------------------------------
rule vsearch_fastq_report:
    conda:
        "envs/tapirs.yaml"
    input:
        "results/02_trimmed/{library}/{sample}.merged.fastq.gz"
    output:
        fqreport = "reports/vsearch/{library}/{sample}_fq_eestats",
        fqreadstats = "reports/vsearch/{library}/{sample}_fq_readstats"
    shell:
        "vsearch \
        --fastq_eestats {input} \
        --output {output.fqreport} ; \
        vsearch \
        --fastq_stats {input} \
        --log {output.fqreadstats} \
        "

#-----------------------------------------------------
# dereplication
#-----------------------------------------------------
rule vsearch_dereplication:
    conda:
        "envs/tapirs.yaml"
    input:
        "results/02_trimmed/{library}/{sample}.merged.fasta"
    output:
        "results/02_trimmed/{library}/{sample}.merged.derep.fasta"
    shell:
        "vsearch \
        --derep_fulllength {input} \
        --sizeout \
        --minuniquesize 3 \
        --output {output} \
        "

#-----------------------------------------------------
# denoise
#-----------------------------------------------------
rule vsearch_denoising:
    conda:
        "envs/tapirs.yaml"
    input:
        "results/02_trimmed/{library}/{sample}.merged.derep.fasta"
    output:
        fasta = "results/03_denoised/{library}/{sample}.fasta",
        biom = "reports/vsearch/{library}/{sample}.denoise.biom"
    #params:
    #    log="reports/denoise/{library}/vsearch.log"
    shell:
        "vsearch \
        --cluster_unoise {input} \
        --centroids {output.fasta} \
        --biomout {output.biom} \
        --minsize 3 \
        --unoise_alpha 0.5 \
        "
        #" --notrunclabels
        # --log {params.log} \

#-----------------------------------------------------
# chimera removal
#-----------------------------------------------------
rule vsearch_dechimerisation: # output needs fixing
    conda:
        "envs/tapirs.yaml"
    input:
        "results/03_denoised/{library}/{sample}.fasta"
    output: # fix
        text = "results/03_denoised/{library}/{sample}_chimera.txt",
        fasta = "results/03_denoised/{library}/nc_{sample}.fasta"
    shell:
        "vsearch \
        --uchime_ref {input} \
        --db data/databases/12S_full/12s_full.fasta \
        --mindiffs 1 \
        --mindiv 0.8 \
        --uchimeout {output.text} \
        --nonchimeras {output.fasta} \
        "

#------------------------------------------------------
# re-replication
#-------------------------------------------------------
rule vsearch_rereplication:
    conda:
        "envs/tapirs.yaml"
    input:
        "results/03_denoised/{library}/{sample}.fasta" # check
    output:
        "results/rereplicated/{library}/{sample}.fasta"
    threads:
        6
    shell:
        "vsearch \
        --rereplicate {input} \
        --output {output} \
        "

#-----------------------------------------------------
# blastn, sequence similarity search
#-----------------------------------------------------
rule blastn:
    #message: "executing blast analsyis of sequences against database {input.database}"
    conda:
        "envs/tapirs.yaml"
    input:
        #db = "nt", #specify in environment.yaml
        query = "results/03_denoised/{library}/nc_{sample}.fasta"
    params:
        db_dir = directory("data/databases/12S_full/12s_full"), # database directory
        descriptions = "50", # return maximum of 50 hits
        outformat = "'6 qseqid stitle sacc staxids pident qcovs evalue bitscore'",
        min_perc_ident = "100", # this needs to be 100%
        min_evalue = "1e-20"
    output:
        "results/blast/{library}/{sample}_blast.out"
    threads:
        6
    shell:
        "blastn \
        -db {params.db_dir} \
        -num_threads {threads} \
        -outfmt {params.outformat} \
        -perc_identity {params.min_perc_ident} \
        -evalue {params.min_evalue} \
        -max_target_seqs {params.descriptions} \
        -query {input.query} \
        -out {output} \
        "

#-----------------------------------------------------
# tax_to_blast, adds taxonomy in a column to blast output
#-----------------------------------------------------
rule tax_to_blast:
    conda:
        "envs/tapirs.yaml"
    input:
        blast_out = "results/blast/{library}/{sample}_blast.out",
        ranked_lineage = "data/databases/new_taxdump/rankedlineage.dmp"
    output:
        blast_taxonomy = "results/blast/{library}/{sample}_tax.tsv"
    shell:
        "python scripts/tax_to_blast.py -i {input.blast_out} -o {output.blast_taxonomy} -lin {input.ranked_lineage}"

#-----------------------------------------------------
# MLCA, majority lowest common ancestor
#-----------------------------------------------------
rule mlca:
    input:
        "results/blast/{library}/{sample}_tax.tsv"
    output:
        "results/mlca/{library}/{sample}_lca.tsv"
    params:
        bitscore = "10", # -b blast hit bitscore upper threshold
        identity = "100", # -id percent identity
        coverage = "60", # -cov percentage coverage
        majority = "100", # -m majority percent, 100 is all hits share taxonomy
        hits = "1" # -hits minimum number of hits, default = 2, 1 is true LCA just takes top hit
    shell:
        "python \
        scripts/mlca.py \
        -i {input} \
        -o {output} \
        -b {params.bitscore} \
        -id {params.identity} \
        -cov {params.coverage} \
        -m {params.majority} \
        -hits {params.hits} \
        "

#-----------------------------------------------------
# Kraken, kmer based taxonomic id
#-----------------------------------------------------
rule kraken2:
    conda:
        "envs/tapirs.yaml"
    input:
        "results/rereplicated/{library}/{sample}.fasta"
    output:
        kraken_outputs = "results/kraken/outputs/{library}.{sample}.tsv",
        kraken_reports = "results/kraken/reports/{library}/{sample}.txt"
    threads:
        6
    params:
        confidence = "0.0",
        kraken_db = directory("data/databases/kraken/kraken2_db") # This is specified but not called in the shell command - Mike
    shell:
        "kraken2 \
        --db data/databases/kraken2_db/ {input} \
        --use-names \
        --memory-mapping \
        --threads {threads} \
        --confidence {params.confidence} \
        --output {output.kraken_outputs} \
        --report {output.kraken_reports} \
        "

# could use --report-zero-counts if against small database
    # will add this t the config file - Mike

#-----------------------------------------------------
# Kraken output to BIOM format
#-----------------------------------------------------
rule kraken_to_biom:
    conda:
        "envs/tapirs.yaml"
    output:
        "results/kraken/{my_experiment}.biom" #my_experiment=config["my_experiment"])
    params:
        input = expand("results/kraken/reports/{library}.{sample}.txt", library=library, sample=sample)
    shell:
        "kraken-biom \
        {params.input} \
        --max F \
        -o {output} \
        "

#---------------------------------------------------
# Biom convert, BIOM to tsv
#---------------------------------------------------
rule biom_convert:
    conda:
        "envs/tapirs.yaml"
    input:
        expand("results/kraken/{my_experiment}.biom", my_experiment=my_experiment)
    output:
        expand("results/kraken/{my_experiment}.tsv", my_experiment=my_experiment)
    threads:
        6
    shell:
        "biom convert -i {input} -o {output} --to-tsv --header-key taxonomy"

#---------------------------------------------------
# sintax, kmer similarity taxonomic ID
#---------------------------------------------------
rule sintax:
    input:
        database = "data/databases/sintax/12s.fas",
        query = "results/rereplicated/{library}/{sample}.fasta"
    output:
        "results/sintax/{library}/{sample}_reads.sintax"
    params:
        cutoff = "0.8"
    shell:
        "vsearch -sintax \
        {input.query} \
        -db {input.database} \
        -tabbedout {output} \
        -strand both \
        -sintax_cutoff {params.cutoff} \
        "

#-------------------------------------------------
# biom taxonomy transformation
#-------------------------------------------------

# rule transform_biomtsv:
#     input:
#         "results/kraken/{my_experiment}.tsv"
#     output:
#         "results/kraken/{my_experiment}.trans.tsv"
#     run:
#         "

#-----------------------------------------------------
# Krona, interactive html graphics of taxonomic diversity
#-----------------------------------------------------
# Multiple files can be given separated by comma.

# requires running : $ ktUpdateTaxonomy
rule kraken_to_krona: # see here: https://github.com/marbl/Krona/issues/117
    conda:
        "envs/tapirs.yaml"
    input:
        "results/kraken/reports/{library}/{sample}.txt"
    output:
        "reports/krona/kraken/{library}/{sample}.html"
    shell:
        "ktImportTaxonomy -m 3 -t 5 {input} -o {output}"

rule mlca_to_krona:
    conda:
        "envs/tapirs.yaml"
    input:
        "results/mlca/{library}/{sample}_lca.tsv"
    output:
        "reports/krona/mlca/{library}/{sample}.html"
    shell:
        "ktImportText -m 2 -t 3 {input} -o {output}"

##this above needs the first column cutting down to just read count - Mike


#-----------------------------------------------------
# convert to tsv format for importing with kronatext
# takes 4th column of input (ie taxonomy passing SINTAX cutoff)
# and exports each unique taxonomy with a count to new tsv
# then rule passes tsv to krona for html plots
#-----------------------------------------------------
rule sintax_to_kronatext:
    input:
        "results/sintax/{library}/{sample}_reads.sintax"
    output:
        "results/sintax/{library}/{sample}_sintax_taxcount.tsv"
    shell:
        "awk '{print $4}' {input} | sort | uniq -c >> {output}"

rule sintaxtext_to_krona: # importing with kronatext to krona
    input:
        "results/sintax/{library}/{sample}_sintax_taxcount.tsv"
    output:
        "reports/krona/{library}/{sample}.sintax.html"
    shell:
        "ktImportText {input} -o {output}"

## DO NOT LOSE THIS COMMAND!!!!
## python /home/mike/anaconda3/pkgs/basta-1.3-py27_1/bin/basta2krona.py
# Desktop/tapirs/results/LCA/testlib/BLEL01.basta_LCA.out Desktop/kronatest.html

#-----------------------------------------------------
# Snakemake report
#-----------------------------------------------------
rule snakemake_report:
    output:
        expand("reports/{my_experiment}_smk-report.html", my_experiment=my_experiment)
    shell:
        "snakemake --report {output}"

#-----------------------------------------------------
# Archive conda environment
#-----------------------------------------------------

rule conda_env:
    conda:
        "envs/{conda_envs}"
    output:
        "reports/archived_envs/{conda_envs}"
    shell:
        "conda env export --file {output}"

##################################################################################################

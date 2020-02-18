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



library="N1"
sample,= glob_wildcards("data/01_demultiplexed/N1/{sample}.R1.fastq.gz")
#sample="BLEL01" #this is temporary while we sort the output of demultiplexing to be zipped
R=["R1", "R2"]
#, = glob_wildcards("data/01_demultiplexed/{library}/")
## check libraries and library are named OK throughout

conda_envs=["tapirs.yaml", "basta_LCA.yaml"]
#sample,= glob_wildcards("data/01_demultiplexed/{library}/{sample}.R1.fastq")

#-----------------------------------------------------
# target rule, specify outputs
#-----------------------------------------------------
rule all:
    input:
        expand("results/02_trimmed/{library}/{sample}.{R}.fastq.gz", library=library, sample=sample, R=R),
        expand("results/03_denoised/{library}/{sample}.fasta", library=library, sample=sample, R=R),
        expand("results/blast/{library}/{sample}_blast.out", library=library, sample=sample),
        #expand("results/LCA/{library}/{sample}.basta_LCA.out", library=library, sample=sample),
        #expand("results/blast/{library}/{sample}_blast.taxed.out", library=library, sample=sample),
        ##expand("results/simpleLCA/{library}/{sample}.lca", library=library, sample=sample),
        #expand("results/LCA/{library}/{sample}.basta_LCA.out.biom", library=library, sample=sample),
		#expand("results/basta/{sample}.basta_LCA.out", library=library, sample=sample),
        expand("results/mlca/{library}/{sample}_lca.tsv", library=library, sample=sample),
        # reports ----------------------------------------------
        expand("reports/fastp/{library}/{sample}.json", library=library, sample=sample),
        expand("reports/fastp/{library}/{sample}.html", library=library, sample=sample),
        expand("reports/vsearch/{library}/{sample}.denoise.biom", library=library, sample=sample),
        expand("reports/vsearch/{library}/{sample}_fq_eestats", library=library, sample=sample),
        expand("reports/vsearch/{library}/{sample}_fq_readstats", library=library, sample=sample),
        #expand("reports/krona/{library}/{sample}.basta_to_krona.html", library=library, sample=sample),
        expand("reports/archived_envs/{conda_envs}", conda_envs=conda_envs),
        #"reports/{my_experiment}_smk-report.html", my_experiment = (config["my_experiment"])
    #    expand("results/LCA/{library}/{sample}.basta_LCA.out.biom", library=library, sample=sample),
        #    expand("results/LCA/{library}/{sample}.basta_LCA.out.tsv", library=library, sample=sample)

#-----------------------------------------------------
# include rule files
#-----------------------------------------------------
# include: "rules/reports.smk",
# include: "rules/kraken.smk",
# include: "rules/blast.smk",
# include: "rules/qc.smk"

#-----------------------------------------------------
# fastp, control for sequence quality and pair reads
#-----------------------------------------------------
rule fastp_trim_and_merge:
    message: "Beginning fastp QC of raw data"
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
        merged="results/02_trimmed/{library}/{sample}.merged.fastq.gz",
        json = "reports/fastp/{library}/{sample}.json",
        html = "reports/fastp/{library}/{sample}.html"
    shell:
        "fastp \
        -i {input.read1} \
        -I {input.read2} \
        -o {output.out1} \
        -O {output.out2} \
        -j {output.json} \
        -h {output.html} \
        --qualified_quality_phred 30 \
        --length_required 90 \
        --unpaired1 {output.out_unpaired1} \
        --unpaired2 {output.out_unpaired2} \
        --failed_out {output.out_failed} \
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
# vsearch, convert files from fastq to fasta
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
# vsearch fastq report
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
# denoise, remove sequence errors
#-----------------------------------------------------
rule vsearch_denoising:
    conda:
        "envs/tapirs.yaml"
    input:
        "results/02_trimmed/{library}/{sample}.merged.derep.fasta"
    output:
        fasta="results/03_denoised/{library}/{sample}.fasta",
        biom="reports/vsearch/{library}/{sample}.denoise.biom"
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
# chimera removal, vsearch
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
        db_dir=directory("data/databases/12S_full/12s_full"), # database directory
        descriptions="50", # return maximum of 50 hits
        outformat="'6 qseqid stitle sacc staxids pident qcovs evalue bitscore'",
        min_perc_ident="100", # this needs to be 100%
        min_evalue="1e-20"
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
        -out {output}"

#-----------------------------------------------------
# LCA, Last Comomon Ancestor analysis of blast using BASTA
#-----------------------------------------------------
# rule basta_LCA:
#     conda:
#         "envs/basta_LCA.yaml"
#     input:
#         "results/blast/{library}/{sample}_blast.out" #fix this
#         # file of blast tabular -outfmt 6 from above
#     params:
#         nhits="50", # -n max number of  hits to consider for classification (default=0=all)
#         minhits="3", # -m must have at least 3 hits, else ignored (default=3)
#         evalue="1e-20", # -e min e-value of hit (default=0.00001)
#         length="90", # -l match must be at least 90bp (default=100)
#         minident="95", # -i minimum identity of hit (default=80)
#         maj_percent="90", # -p 90 = taxonomy shared by 9/10 hits, (default=100 = shared by all)
#         dir="/media/mike/mikesdrive/" # -d directory of database files (default: $HOME/.basta/taxonomy)
#     output: # check library/sample syntax
#         "results/LCA/{library}/{sample}.basta_LCA.out"
#     shell:
#         "basta sequence {input} {output} gb \
#         -p {params.maj_percent} \
#         -m {params.minhits} \
#         -l {params.length} \
#         -i {params.minident} \
#         -n {params.nhits}"
# #        "./bin/basta multiple INPUT_DIRECTORY OUTPUT_FILE MAPPING_FILE_TYPE"
#
#-----------------------------------------------------
# tax_to_blast, adds taxonomy in a column to blast output
#-----------------------------------------------------
rule tax_to_blast:
    conda:
        "envs/tapirs.yaml"
    input:
        blast_out="results/blast/{library}/{sample}_blast.out",
        ranked_lineage="data/databases/new_taxdump/rankedlineage.dmp"
    output:
        blast_taxonomy="results/blast/{library}/{sample}_tax.tsv"
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
        -hits {params.hits}"

#-----------------------------------------------------
# Kraken, kmer based taxonomic id
#-----------------------------------------------------
rule kraken2:
  input:
    seqs = "results/rereplicated/{library}/{sample}.fasta",
    kraken_db = directory("data/databases/kraken/kraken2_db")
  output:
    kraken_outputs = "/results/kraken/outputs/{sample}.tsv",
    kraken_reports = "results/kraken/report/{sample}.txt"
params:
    threads="6",
    confidence="0.0"
  shell:
    "kraken2 --db fish_db {input.seqs} \
    --use-names \
    --memory-mapping \
    --threads {params.threads} \
    --confidence {params.confidence} \
    --output {output.kraken_outputs} \
    --report {output.kraken_reports} "
# could use --report-zero-counts if against small database

#-----------------------------------------------------
# Kraken output to BIOM format
#-----------------------------------------------------
rule kraken_to_biom:
    input:
        directory("results/kraken/report/")
    output:
        "results/kraken/{my_experiment}.biom", my_experiment=config["my_experiment"]
    shell:
        "kraken-biom {input} -max F -o {output}"

#-----------------------------------------------------
# Krona, interactive html graphics of taxonomic diversity
#-----------------------------------------------------
# Multiple files can be given separated by comma.

rule kraken_to_krona: # see here: https://github.com/marbl/Krona/issues/117
    input:
        kraken_output = "/results/kraken/outputs/{sample}.tsv"
    output:
        "reports/krona/kraken/{sample}.html",
    script:
        "scripts/krona/ImportTaxonomy.pl -q 2 -t 3 {input.kraken_output} -o {output}"

rule mlca_to_krona:
    input:
        "results/mlca/{library}/{sample}_lca.tsv"
    output:
        "reports/krona/krona_mlca.html"
    script: # check grammar
        "scripts/krona/ImportText.pl {input} -o {output}"

## DO NOT LOSE THIS COMMAND!!!!
## python /home/mike/anaconda3/pkgs/basta-1.3-py27_1/bin/basta2krona.py
# Desktop/tapirs/results/LCA/testlib/BLEL01.basta_LCA.out Desktop/kronatest.html

#-----------------------------------------------------
# Snakemake report
#-----------------------------------------------------
rule snakemake_report:
    output:
        "reports/{my_experiment}_smk-report.html", my_experiment = (config["my_experiment"])
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

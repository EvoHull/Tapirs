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
        expand("results/LCA/{library}/{sample}.basta_LCA.out", library=library, sample=sample),
        #expand("results/blast/{library}/{sample}_blast.taxed.out", library=library, sample=sample),
        ##expand("results/simpleLCA/{library}/{sample}.lca", library=library, sample=sample),
        #expand("results/LCA/{library}/{sample}.basta_LCA.out.biom", library=library, sample=sample),
		#expand("results/basta/{sample}.basta_LCA.out", library=library, sample=sample),
# reports ----------------------------------------------
        expand("reports/fastp/{library}/{sample}.json", library=library, sample=sample),
        expand("reports/fastp/{library}/{sample}.html", library=library, sample=sample),
        expand("reports/vsearch/{library}/{sample}.denoise.biom", library=library, sample=sample),
        expand("reports/vsearch/{library}/{sample}_fq_eestats", library=library, sample=sample),
        expand("reports/vsearch/{library}/{sample}_fq_readstats", library=library, sample=sample),
        expand("reports/krona/{library}/{sample}.basta_to_krona.html", library=library, sample=sample),
        expand("reports/archived_envs/{conda_envs}", conda_envs=conda_envs),
    #    expand("results/LCA/{library}/{sample}.basta_LCA.out.biom", library=library, sample=sample),
        #    expand("results/LCA/{library}/{sample}.basta_LCA.out.tsv", library=library, sample=sample)

#-----------------------------------------------------
# include rule files
#-----------------------------------------------------
# include: "rules/reports.smk",
# include: "rules/kraken.smk",
# include: "rules/blast.smk",
# include: "rules/qc.smk"


# #-----------------------------------------------------
# # gzip demultiplexed files, seqkit
# # should modify demultiplex.py to do this
# #-----------------------------------------------------
# rule gzip:
#     input:
#         "data/01_demultiplexed/{library}/{sample}.{R}.fastq"
#     output:
#         "data/1_demultiplexed/{library}/{sample}.{R}.fastq.gz"
#     shell:
#         "gzip {input} > {output}"


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
        --uchime3_denovo {input} \
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
        outformat="'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'",
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

# database is going to cause problems and needs a symbolic path in config

#-----------------------------------------------------
# LCA, Last Comomon Ancestor analysis of blast using BASTA
#-----------------------------------------------------
rule basta_LCA:
    conda:
        "envs/basta_LCA.yaml"
    input:
        "results/blast/{library}/{sample}_blast.out" #fix this
        # file of blast tabular -outfmt 6 from above
    params:
        nhits="50", # -n max number of  hits to consider for classification (default=0=all)
        minhits="3", # -m must have at least 3 hits, else ignored (default=3)
        evalue="1e-20", # -e min e-value of hit (default=0.00001)
        length="90", # -l match must be at least 90bp (default=100)
        minident="95", # -i minimum identity of hit (default=80)
        maj_percent="90", # -p 90 = taxonomy shared by 9/10 hits, (default=100 = shared by all)
        dir="/media/mike/mikesdrive/" # -d directory of database files (default: $HOME/.basta/taxonomy)
    output: # check library/sample syntax
        "results/LCA/{library}/{sample}.basta_LCA.out"
    shell:
        "basta sequence {input} {output} gb \
        -p {params.maj_percent} \
        -m {params.minhits} \
        -l {params.length} \
        -i {params.minident} \
        -n {params.nhits}"
#        "./bin/basta multiple INPUT_DIRECTORY OUTPUT_FILE MAPPING_FILE_TYPE"

#-----------------------------------------------------
# Simple-LCA
#-----------------------------------------------------
# rule simpleLCA_adding_taxid:
#     conda:
#         "envs/tapirs.yaml"
#     input:
#         "results/blast/{library}/{sample}_blast.out"
#     output:
#         "results/blast/{library}/{sample}_blast.taxed.out"
#     threads:
#         28
#     shell:
#         "scripts/Simple-LCA-master/add_taxonomy.py -i {input} -t rankedlineage.dmp -m merged.dmp -o {output}"
#
# rule simpleLCA:
#     conda:
#         "envs/tapirs.yaml"
#     input:
#         "results/blast/{library}/{sample}_blast.taxed.out"
#     params:
#         bitscore = "8", # '-b', '--bitscore', 'bitscore top percentage threshold',
#         id = "80", # -id identity threshold, required
#         coverage = "80", # -cov coverage threshold', required=True
#         tophit = "yes", # -t 'Check the best hit first, if it is above the given threshold the tophit will become the output', required=False, choices=['no', 'yes'], default='no')
#         tid = "99", # 'identity treshold for the tophit', required=False, default='100'
#         tcov = "100", # query coverage threshold for the tophit', required=False,  default='100'
#         fh = "environmental", # filter out hits that contain unwanted taxonomy'
#         flh = "unknown", # -flh', filter lca hits', dest='filterLcaHits', help='ignore this in the lca determination'
#     output:
#         "results/simpleLCA/{library}/{sample}.lca"
#     shell:
#         "scripts/Simple-LCA-master/lca.py \
#         -i {input} \
#         -o {output} \
#         -b {params.bitscore} \
#         -id {params.id} \
#         -cov {params.coverage} \
#         -t {params.tophit} \
#         -tid {params.tid} \
#         -tcov {params.tcov} \
#         -fh {params.fh} \
#         -flh {params.flh}"
#         # "scripts/Simple-LCA-master/lca.py -i {input} -o {output} -b 8 -id 80 -cov 80 -t yes -tid 99 -tcov 100 -fh 'environmental' -flh 'unknown'"

#-----------------------------------------------------
# BASTA to BIOM,
#-----------------------------------------------------
# BASTA output tsv converted to BIOM, uses BIOM-convert

# rule basta_BIOM:
#     conda:
#         "envs/tapirs.yaml"
#     input:
#         "results/LCA/{library}/{sample}.basta_LCA.out"
#     params:
#         json="json",
#         hdf5="hdf5"
#     output:
#         "results/LCA/{library}/{sample}.basta_LCA.out.biom"
#     shell:
#         "biom convert \
#         -i {input} \
#         -o {output} \
#         --table-type='OTU table' \
#         --to-{params.hdf5} \
#         "

#-----------------------------------------------------
# BIOM to tsv GRAHAM TO CHECK
#-----------------------------------------------------
# rule BIOM_tsv:
#     conda:
#         "envs/tapirs.yaml"
#     input:
#         "results/LCA/{library}/{sample}.basta_LCA.out.biom"
#     output:
#         "results/LCA/{library}/{sample}.basta_LCA.out.tsv"
#     shell:
#         "biom convert -i {input} -o {output} --table-type='OTU table' --to-{params.hdf5}"


# biom convert -i table.txt -o table.from_txt_json.biom --table-type="OTU table" --to-json
# biom convert -i table.txt -o table.from_txt_hdf5.biom --table-type="OTU table" --to-hdf5
# OUTPUT: workflow should export data for downstream analysis. This is BIOM format written by metaBEAT, and also csv I guess.

#-----------------------------------------------------
# Krona
#-----------------------------------------------------
# basta2krona.py
# This creates a krona plot (html file) for each sample that can be opened in a browser from a basta annotation output file(s).
# Multiple files can be given separated by comma.

rule krona_LCA_plot:
    conda:
        "envs/tapirs.yaml"
    input:
        "results/LCA/{library}/{sample}.basta_LCA.out"
    output:
        "reports/krona/{library}/{sample}.basta_to_krona.html"
    params:
        path="/home/mike/anaconda3/pkgs/basta-1.3-py27_1/bin" # this is going to need sorting out in the config file
    shell:
        "python \
        {params.path}/basta2krona.py \
        {input} \
        {output}"



## DO NOT LOSE THIS COMMAND!!!!
## python /home/mike/anaconda3/pkgs/basta-1.3-py27_1/bin/basta2krona.py
# Desktop/tapirs/results/LCA/testlib/BLEL01.basta_LCA.out Desktop/kronatest.html


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

# ==================================================
# TAPIRS REPORT GENERATION
# ==================================================
# A workflow reporting on QC and taxonomic assignment
# fastp reports are written from the qc.smk rule

configfile: "config.yaml"
# report: "reports/snakemake-report.html"

# --------------------------------------------------
# Snakemake, report
# --------------------------------------------------

# rule snakemake_report:  # this only works when run after the workflow has completed, so unsure how to work it into the snakemake
#     output:
#         "reports/snakemake-report.html"
#     shell:
#         "snakemake --report {output}"

# --------------------------------------------------
# Snakemake, plot DAG figure of workflow
# --------------------------------------------------

rule plot_workflow_DAG:
    output:
        "reports/dag_rulegraph.png"
    shell:
        "snakemake --rulegraph | dot -Tpng > {output}"

# --------------------------------------------------
# Conda, archive environment
# --------------------------------------------------

rule conda_env:
    output:
        "reports/archived_envs/tapirs.yaml"
    shell:
        "conda env export --file {output}"

# --------------------------------------------------
# Seqkit, write report on 02_trimmed files
# for each library make a report on all files
# --------------------------------------------------

# rule seqkit_stats_trimmedfiles:
#     input:
#         "results/02_trimmed/{SAMPLES}"
#     threads:
#         4  # -j
#     output:
#         tsv = "reports/seqkit/{SAMPLES}.trimmed.seqkit-stats.tsv",
#         md = "reports/seqkit/{SAMPLES}.trimmed.seqkit-stats.md",
#     shell:
#         """
#         seqkit stats {input}*.fastq -b -e -T -j 4 -o {output.tsv} ;
#         csvtk csv2md {output.tsv} -t -o {output.md}
#         """

# -----------------------------------------------------
# vsearch, readstats report on 03_merged concatenated fastq
# -----------------------------------------------------

rule vsearch_fastq_readstats:
    conda:
        "../envs/environment.yaml"
    input:
        expand("results/03_merged/{SAMPLES}.concat.fastq", SAMPLES=SAMPLES)
    output:
        fqreadstats = "reports/vsearch/{SAMPLES}.concat.fq_readstats"
    shell:
        """
        vsearch --fastq_stats {input} --log {output.fqreadstats}
        """

# -----------------------------------------------------
# vsearch, eestats report on 03_merged concatenated fastq
# -----------------------------------------------------

rule vsearch_fastq_eestats:
    conda:
        "../envs/environment.yaml"
    input:
        expand("results/03_merged/{SAMPLES}.concat.fastq", SAMPLES=SAMPLES)
    output:
        fqreport = "reports/vsearch/{SAMPLES}.concat.fq_eestats",
    shell:
        """
        vsearch --fastq_eestats {input} --output {output.fqreport};
        """

#---------------------------------------------------
# Seqkit, report on 03_merged concatenated fasta
#---------------------------------------------------

# rule seqkit_stats_mergedfiles:
#     input:
#         expand("results/03_merged/{library}/{sample}.concat.fasta", sample=SAMPLES, library=LIBRARIES)
#     threads:
#         12
#     output:
#         tsv = "reports/seqkit/{library}.concat.seqkit-stats.tsv",
#         md = "reports/seqkit/{library}.concat.seqkit-stats.md",
#     shell:
#         """
#         seqkit stats {input} -b -e -T -j {threads} -o {output.tsv} ;
#         csvtk csv2md {output.tsv} -t -o {output.md} ;
#         """

#---------------------------------------------------
# MultiQC, aggregate QC reports as html report
#-----------------------------------------------------

# rule multiqc:
#     conda:
#         "envs/environment.yaml"
#     input:
#         "reports/fastp/{library}/"
#     output:
#         filename = "{library}.multiqc.html",  # report filename
#         outdir = directory("reports/multiqc/")
#     params:
#         overwrite = "-f",  # overwrite previous multiqc output
#         zip = "-z",  # compress the multiqc data dir
#         quiet = "-q", # only log errors
#         dirnames = "-dd 1",  # prepend library dir name to sample names
#     shell:
#         "multiqc {input} \
#           -n {output.filename} \
#           -o {output.outdir} \
#           {params.dirnames} \
#           {params.overwrite} \
#           {params.zip} \
#           {params.quiet} \
#           "

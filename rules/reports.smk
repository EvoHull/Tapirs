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

rule snakemake_report:
    output:
        "reports/snakemake-report.html"
    shell:
        "snakemake --report {output}"

# --------------------------------------------------
# Snakemake, report from docs
# --------------------------------------------------

rule report:
    input:
        "reports/rulegraph_dag.png"
    output:
        "snakemake-report.html"
    run:
        from snakemake.utils import report
        # with open(input[0]) as vcf:
        #     n_calls = sum(1 for l in vcf if not l.startswith("#"))

        report(output[0])

# --------------------------------------------------
# Snakemake, plot DAG figure of workflow
# --------------------------------------------------

rule plot_workflow_DAG:
    output:
        "reports/rulegraph_dag.png"
    shell:
        "snakemake --rulegraph | dot -Tpng > {output}"

# --------------------------------------------------
# Conda, archive environment
# --------------------------------------------------

rule conda_env:
    conda:
        "../envs/{conda_envs}"
    output:
        "reports/archived_envs/{conda_envs}"
    shell:
        "conda env export --file {output}"

# --------------------------------------------------
# Seqkit, write report on 02_trimmed files
# for each library make a report on all files 
# --------------------------------------------------

rule seqkit_stats_trimmedfiles:
    input:
        expand("results/02_trimmed/{sample.library}/",
            sample=sample.reset_index().itertuples()),
    threads:
        12
    output:
        tsv = "reports/seqkit/{library}.trimmed.seqkit-stats.tsv",
        md = "reports/seqkit/{library}.trimmed.seqkit-stats.md",
        # histogram = report("reports/seqkit/{library}.av-length-histogram", category="QC")
        # report("fig1.svg", caption="report/fig1.rst", category="Step 1")
    shell:
        """
        seqkit stats {input}/* -b -e -T -j {threads} -o {output.tsv} ;
        csvtk csv2md {output.tsv} -t -o {output.md} ;
        """ 
#             csvtk -t plot hist {output.tsv} -f 6 -o {output.histogram}

# -----------------------------------------------------
# vsearch, report on 03_merged concatenated fastq
# -----------------------------------------------------

rule vsearch_fastq_report:
    conda:
        "../envs/environment.yaml"
    input:
        # "results/03_merged/{library}/{sample}.concat.fastq"
        expand("results/03_merged/{sample.library}/{sample.sample}.concat.fastq", 
            sample=sample.reset_index().itertuples()),
    output:
        fqreport = "reports/vsearch/{library}/{sample}.concat.fq_eestats",
        fqreadstats = "reports/vsearch/{library}/{sample}.concat.fq_readstats"
    shell:
        "vsearch --fastq_eestats {input} --output {output.fqreport}"
        "vsearch --fastq_stats {input} --log {output.fqreadstats}"

#---------------------------------------------------
# Seqkit, report on 03_merged concatenated fasta
#---------------------------------------------------

rule seqkit_stats_mergedfiles:
    input:
        # "results/03_merged/{library}/{sample}.concat.fasta"
        expand("results/03_merged/{sample.library}/{sample.sample}.concat.fasta", sample=sample.reset_index().itertuples()),
    threads:
        12
    output:
        tsv = "reports/seqkit/{library}.concat.seqkit-stats.tsv",
        md = "reports/seqkit/{library}.concat.seqkit-stats.md",
    shell:
        """
        seqkit stats {input}/* -b -e -T -j {threads} -o {output.tsv} ;
        csvtk csv2md {output.tsv} -t -o {output.md} ;
        """ 

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

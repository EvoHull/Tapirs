# ==================================================
# TAPIRS REPORT GENERATION
# ==================================================
# A workflow reporting on QC and taxonomic assignment
# Some fastp reports are written from the qc.smk rule

configfile: "config.yaml"
report: "reports/snakemake-report.html"
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
# Seqkit, write simple report on 02_trimmed files
# --------------------------------------------------

rule seqkit_stats_allfiles:
        input:
            "results/02_trimmed/{library}"
        threads:
            12
        output:
            tsv = "reports/seqkit/{library}_seqkit-stats.tsv",
            md = "reports/seqkit/{library}_seqkit-stats.md",
            histogram = report("reports/seqkit/{library}_av-length.jpg", category="QC")
            # report("fig1.svg", caption="report/fig1.rst", category="Step 1")
        shell:
            """
            seqkit stats {input}/* -b -e -T -j {threads} -o {output.tsv} ;
            csvtk csv2md {output.tsv} -t -o {output.md} ;
            csvtk -t plot hist {output.tsv} -f 6 -o {output.histogram}
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
#         "reports/multiqc/{library}.multiqc.html"
#     params:
#         outdir = directory("reports/multiqc/"),  # location for report
#         filename = "{library}.multiqc.html",  # report filename
#         overwrite = "-f",  # overwrite previous multiqc output
#         zip = "-z",  # compress the multiqc data dir
#         quiet = "-q", # only log errors
#         dirnames = "-dd 1",  # prepend library dir name to sample names
#     shell:
#         "multiqc {input} -n {params.filename} {params.overwrite} {params.dirnames} {params.zip} {params.quiet} -o {params.outdir}"

# --------------------------------------------------
# Vegan
# --------------------------------------------------
# add vegan graphical and statistical outputs

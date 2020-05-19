# ==================================================
# TAPIRS REPORT GENERATION
# ==================================================
# A workflow reporting on QC and taxonomic assignment
# Some fastp reports are written from the qc.smk rule

configfile: "config.yaml"
report: "../reports/tapirs.rst"


# --------------------------------------------------
# Snakemake, report
# --------------------------------------------------

rule snakemake_report:
    output:
        "reports/snakemake-report.html"
    shell:
        "snakemake --report {output}"


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
# Seqkit, write simple report on fasta files
# --------------------------------------------------

rule seqkit_stats:
        input:
            "results/02_trimmed/{library}/{sample}.merged.fasta"
        output:
            "reports/seqkit/{library}/{sample}.fastastats.md"
        shell:
            "seqkit stats {input} | csvtk csv2md -t -o {output}"


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

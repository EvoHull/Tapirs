# ==================================================
# TAPIRS REPORT GENERATION
# ==================================================
# A workflow reporting on QC and taxonomic assignment

configfile: "config.yaml"
report: "../reports/tapirs.rst"

# --------------------------------------------------
# seqkit, write simple report on fasta files
# --------------------------------------------------
rule seqkit_stats:
        input:
            "results/02_trimmed/{library}/{sample}.merged.fasta"
        output:
            "reports/seqkit/{library}/{sample}.fastastats.md"
        shell:
            "seqkit stats {input} | csvtk csv2md -t -o {output}"


# --------------------------------------------------
# Vegan
# --------------------------------------------------

# --------------------------------------------------
# Snakemake report
# --------------------------------------------------
rule snakemake_report:
    output:
        expand("reports/{my_experiment}_smk-report.html", my_experiment=config["my_experiment"])
    shell:
        "snakemake --report {output}"


# --------------------------------------------------
# Archive conda environment
# --------------------------------------------------

rule conda_env:
    conda:
        "../envs/{conda_envs}"
    output:
        "reports/archived_envs/{conda_envs}"
    shell:
        "conda env export --file {output}"


#---------------------------------------------------
# MultiQC, aggregate QC reports as html report
#-----------------------------------------------------
rule multiqc:
    conda:
        "../envs/environment.yaml"
    input:
        "reports/fastp/{library}/"
    output:
        "reports/multiqc/{library}.multiqc.html"
    params:
        outdir = directory("reports/multiqc/"),  # location for report
        filename = "{library}.multiqc.html",  # report filename
        overwrite = "-f",  # overwrite previous multiqc output
        zip = "-z",  # compress the multiqc data dir
        quiet = "-q", # only log errors
        dirnames = "-dd 1",  # prepend library dir name to sample names
    shell:
        "multiqc {input} \
        -n {params.filename} \
        {params.overwrite} \
        {params.dirnames} \
        {params.zip} \
        {params.quiet} \
        -o {params.outdir}"

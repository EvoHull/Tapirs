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
# MultiQC
#-----------------------------------------------------
rule multiqc:
    conda:
        "../envs/tapirs.yaml"
    input:
        "reports/fastp/{library}/"
    output:
        "reports/multiqc/{library}.multiqc.html"
    params:
        out = "reports/multiqc/",
        n = "{library}.multiqc.html"
    shell:
        "multiqc {input} -n {params.n} -o {params.out}"

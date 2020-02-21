#-----------------------------------------------------
# Tapirs-report
# -------------
# A metabarcoding workflow to report on QC and
# taxonomic assignment
#-----------------------------------------------------

configfile: "config.yaml"

report: "../reports/tapirs.rst"

#-----------------------------------------------------
rule all_reports:
    input:
        expand("data/{sample}.R[1,2].fastq.gz", sample=SAMPLES)

#-----------------------------------------------------
# seqkit to write simple report on fasta files
#-----------------------------------------------------
rule seqkit_stats:
        input:
            "results/02_trimmed/{library}/{sample}.merged.fasta"
        output:
            "reports/seqkit/{library}/{sample}.fastastats.md"
        shell:
            "seqkit stats {input} | csvtk csv2md -t -o {output}"

#-----------------------------------------------------
# Vegan
#-----------------------------------------------------

#-----------------------------------------------------
# Snakemake report
#-----------------------------------------------------
rule snakemake_report:
    output:
        expand("reports/{my_experiment}_smk-report.html", my_experiment=config["my_experiment"])
    shell:
        "snakemake --report {output}"

#-----------------------------------------------------
# Archive conda environment
#-----------------------------------------------------

rule conda_env:
    conda:
        "../envs/{conda_envs}"
    output:
        "reports/archived_envs/{conda_envs}"
    shell:
        "conda env export --file {output}"

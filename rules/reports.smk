#-----------------------------------------------------
# Tapirs-report
# -------------
# A metabarcoding workflow to report on QC and
# taxonomic assignment
#-----------------------------------------------------

report: "../reports/tapirs.rst"

# get the sequence files into a list
SAMPLES, = glob_wildcards("{sample}.R[1,2].fastq.gz")

#-----------------------------------------------------
rule all_reports:
    input:
        expand("data/{sample}.R[1,2].fastq.gz", sample=SAMPLES)


#-----------------------------------------------------
# seqkit to write simple report on fasta files
#-----------------------------------------------------
rule seqkit_stats:
        input:
            "results/seqkit/{sample}.fasta.gz"
        output:
            "reports/{sample}/seqkit_fastastats.md"
        shell:
            "seqkit stats {input} | csvtk csv2md -t -o {output}"

#-----------------------------------------------------
# tidy fastp report
#-----------------------------------------------------
rule fastpreport_mover: # moving reports from fastp
    input:
        "fastp.*" # check this syntax
    output:
        directory("reports/")
    shell:
        "mv {input} {output}"

#-----------------------------------------------------
# Vegan
#-----------------------------------------------------


#-----------------------------------------------------
# Krona
#-----------------------------------------------------
# basta2krona.py
# This creates a krona plot (html file) that can be opened in your browser from a basta annotation output file(s). Multiple files can be given separated by comma.

rule basta_to_krona:
    input:
        "results/basta/basta_LCA.out"
    output:
        "reports/krona/basta_to_krona.html"
    shell:
        "./scripts/basta2krona {input} {output}"

#-----------------------------------------------------
# Snakemake report
#-----------------------------------------------------
rule snakemake_report:
    output:
        "reports/snakemake_report.html"
    shell:
        "snakemake --report {output}"

#-----------------------------------------------------
# Archive conda environment
#-----------------------------------------------------
rule conda_env:
    shell:
        "conda env export > reports/conda_environment.yml"

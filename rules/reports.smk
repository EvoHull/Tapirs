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
# Generate files suitable for Krona reports
#       only one of these rules is required perhaps?
#-----------------------------------------------------
rule biom_to_krona_tsv:
# convert BIOM output of LCA analysis to krona
# https://github.com/GenomicaMicrob/OTUsamples2krona
    input:
        biom = "{sample}_otu_table.biom"
    output:
        tmp_table = temp("results/blast/{sample}_otu_table.tsv"),
        otu_final = "results/blast/{sample}_otu.tsv"
    shell:
        "biom convert -i {input.biom} -o {output.tmp_table} --to-tsv --header-key taxonomy;\
        cut -f1 --complement {output.tmp_table} > {output.otu_final}" # fixes formatting

rule lca_to_krona: # pseudocode
    input:
        "results/blast/lca/{sample}.tsv"
    output:
        "results/blast/lca/{sample}_krona.tsv"
    scriot:
        "scripts/lca_to_krona.py"

#-----------------------------------------------------
# Generate Krona html files
#-----------------------------------------------------
rule basta_to_krona:
# basta2krona.py work from a basta annotation output file(s)
# Multiple files can be given separated by comma.
    input:
        "results/basta/basta_LCA.out"
    output:
        "reports/krona/basta/{sample}.html"
    script:
        "scripts/basta2krona.py {input} {output}"

rule kraken_to_krona: # see here: https://github.com/marbl/Krona/issues/117
    input:
    output:
        "reports/krona/kraken/{sample}.html"
    script:
        "scripts/krona/ImportTaxonomy.pl -q 2 -t 3 YourKrakenOutputFile -o YourKronaReportFile"

rule blastlca_to_krona:
    input:
        "results/blast/{sample}_otu.tsv" # or lca tsv
    output:
        "reports/krona/blast/{sample}.html"
    shell:

#-----------------------------------------------------
# Snakemake report
#-----------------------------------------------------
rule snakemake_report:
    output:
        "reports/(config["my_experiment"])_smk-report.html"
    shell:
        "snakemake --report {output}"

#-----------------------------------------------------
# Archive conda environment
#-----------------------------------------------------
rule conda_env:
    shell:
        "conda env export > reports/conda_environment.yml"

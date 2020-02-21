#####   KRAKEN   #####
#-----------------------------------------------------
# Kraken, kmer based taxonomic id
#-----------------------------------------------------

configfile: "config.yaml"

rule kraken2:
    conda:
        "../envs/tapirs.yaml"
    input:
        "results/rereplicated/{library}/{sample}.fasta"
    output:
        kraken_outputs = "results/kraken/outputs/{library}/{sample}.tsv",
        kraken_reports = "results/kraken/reports/{library}.{sample}.txt"     #same here - Mike
    threads:
        6
    params:
        confidence = "0.0",
        kraken_db = directory("data/databases/kraken/kraken2_db") # This is specified but not called in the shell command - Mike
    shell:
        "kraken2 \
        --db data/databases/kraken2_db/ {input} \
        --use-names \
        --memory-mapping \
        --threads {threads} \
        --confidence {params.confidence} \
        --output {output.kraken_outputs} \
        --report {output.kraken_reports} \
        "

# could use --report-zero-counts if against small database
    # will add this t the config file - Mike


#-----------------------------------------------------
# Kraken output to BIOM format
#-----------------------------------------------------

rule kraken_to_biom:
    conda:
        "../envs/tapirs.yaml"
    output:
        "results/kraken/{my_experiment}.biom" #my_experiment=config["my_experiment"])
    params:
        input = expand("results/kraken/reports/{library}.{sample}.txt", library=library, sample=sample)
    shell:
        "kraken-biom \
        {params.input} \
        --max F \
        -o {output} \
        "


#---------------------------------------------------
# Biom convert, BIOM to tsv
#---------------------------------------------------

rule biom_convert:
    conda:
        "../envs/tapirs.yaml"
    input:
        expand("results/kraken/{my_experiment}.biom", my_experiment=config["my_experiment"])
    output:
        expand("results/kraken/{my_experiment}.tsv", my_experiment=config["my_experiment"])
    threads:
        6
    shell:
        "biom convert -i {input} -o {output} --to-tsv --header-key taxonomy"

#-------------------------------------------------
# biom taxonomy transformation
#-------------------------------------------------

# rule transform_biomtsv:
#     input:
#         "results/kraken/{my_experiment}.tsv"
#     output:
#         "results/kraken/{my_experiment}.trans.tsv"
#     run:
#         "


#-----------------------------------------------------
# Krona, interactive html graphics of taxonomic diversity
#-----------------------------------------------------

rule kraken_to_krona: # see here: https://github.com/marbl/Krona/issues/117
    conda:
        "../envs/tapirs.yaml"
    input:
        "results/kraken/outputs/{library}/{sample}.tsv"
    output:
        "reports/krona/kraken/{library}/{sample}.html"
    shell:
        "ktImportTaxonomy -q 2 -t 3 {input} -o {output}"

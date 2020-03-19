# ==================================================
# SINTAX ANALYSIS
# ==================================================

configfile: "config.yaml"

# convert to tsv format for importing with kronatext
# takes 4th column of input (ie taxonomy passing SINTAX cutoff)
# and exports each unique taxonomy with a count to new tsv
# then rule passes tsv to krona for html plots

# --------------------------------------------------
# sintax, kmer similarity taxonomic ID
# --------------------------------------------------

rule sintax:
    conda:
        "../envs/environment.yaml"
    input:
        query = "results/rereplicated/{library}/{sample}_rerep.fasta"
    output:
        "results/sintax/{library}/{sample}_reads.sintax"
    params:
        database = config["sintax_db"]
    shell:
        "vsearch -sintax \
        {input.query} \
        -db {params.database} \
        -tabbedout {output} \
        -strand both \
        -sintax_cutoff {config[SINTAX_cutoff]} \
        "


# --------------------------------------------------
# convert sintax output to tsv compatible with krona
# --------------------------------------------------

rule sintax_to_krona_transformation:
    conda:
        "../envs/environment.yaml"
    input:
        "results/sintax/{library}/{sample}_reads.sintax"
    output:
        "results/sintax/{library}/{sample}_sintax_taxcount.tsv"  ### this needs transfromation before itll go into krona - Mike
    shell:
        "awk '{{print $4}}' {input} | \
        sort | \
        uniq -c | \
        sed -e 's/^[ \t]*//' \
        | sed -e 's/ /\t/g' | \
        sed -e 's/:/\t/g' | \
        sed -e 's/,/\t/g' | \
        cut -f 1,3,5,7,9,11,13,15 >> {output}"


# --------------------------------------------------
# Sintax to krona
# --------------------------------------------------

rule sintax_krona_plot:
    conda:
        "../envs/environment.yaml"
    input:
        "results/sintax/{library}/{sample}_sintax_taxcount.tsv"
    output:
        "reports/krona/sintax/{library}/{sample}.sintax.html"
    shell:
        "ktImportText {input} -o {output}"

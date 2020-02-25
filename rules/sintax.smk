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
        "../envs/tapirs.yaml"
    input:
        query = "results/rereplicated/{library}/{sample}.fasta"
    output:
        "results/sintax/{library}/{sample}_reads.sintax"
    params:
        cutoff = "0.8",
        database = config["sintax_db"]

    shell:
        "vsearch -sintax \
        {input.query} \
        -db {params.database} \
        -tabbedout {output} \
        -strand both \
        -sintax_cutoff {params.cutoff} \
        "


# --------------------------------------------------
# convert sintax output to tsv compatible with krona
# --------------------------------------------------

rule sintax_to_kronatext:
    conda:
        "../envs/tapirs.yaml"
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

rule sintaxtext_to_krona:
    conda:
        "../envs/tapirs.yaml"
    input:
        "results/sintax/{library}/{sample}_sintax_taxcount.tsv"
    output:
        "reports/krona/sintax/{library}/{sample}.sintax.html"
    shell:
        "ktImportText {input} -o {output}"

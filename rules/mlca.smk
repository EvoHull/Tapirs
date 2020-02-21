#####   MLCA   #####
#-----------------------------------------------------
# MLCA, majority lowest common ancestor
#-----------------------------------------------------

configfile: "config.yaml"

rule mlca:
    input:
        "results/blast/{library}/{sample}_tax.tsv"
    output:
        "results/mlca/{library}/{sample}_lca.tsv"
    params:
        bitscore = "10", # -b blast hit bitscore upper threshold
        identity = "100", # -id percent identity
        coverage = "60", # -cov percentage coverage
        majority = "100", # -m majority percent, 100 is all hits share taxonomy
        hits = "1" # -hits minimum number of hits, default = 2, 1 is true LCA just takes top hit
    shell:
        "python \
        scripts/mlca.py \
        -i {input} \
        -o {output} \
        -b {params.bitscore} \
        -id {params.identity} \
        -cov {params.coverage} \
        -m {params.majority} \
        -hits {params.hits} \
        "

#-----------------------------------------------------
# Krona, interactive html graphics of taxonomic diversity
#-----------------------------------------------------

rule mlca_to_krona:
    conda:
        "../envs/tapirs.yaml"
    input:
        "results/mlca/{library}/{sample}_lca.tsv"
    output:
        "reports/krona/mlca/{library}/{sample}.html"
    shell:
        "ktImportText {input} -o {output}"

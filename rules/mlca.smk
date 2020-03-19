# ==================================================
# LCA: LOWEST COMMON ANCESTOR ANALYSIS
# ==================================================

configfile: "config.yaml"

# --------------------------------------------------
# MLCA, majority lowest common ancestor
# --------------------------------------------------

rule mlca:
    input:
        "results/blasttax/{library}/{sample}_tax.tsv"
    output:
        "results/mlca/{library}/{sample}_lca.tsv"
    params:
        out = "results/mlca/{library}/{sample}_lca.tsv",
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
# GS - The mlca script needs changing because o fan error. The final species name is output seperated by a tab and not an undderscore or space; is this something you can fix?

#-------------------------------------------------------
# Mlca to tsv
#---------------------------------------------------------

# rule mlca2tsv_transform:
#     input:
#         expand("results/mlca/{sample.library}/{sample.sample}_lca.tsv", sample=sample.reset_index().itertuples())
#     output:
#         directory("reports/mlca/mlcatmp/")
#     params:
#         "reports/mlca/mlcatmp/"
#     shell:
#         "rm -rf {params} \
#         && mkdir {params} \
#         && cp {input} {params}"

rule mlca2tsv:
    # input:
    #     "results/mlca/{library}/{sample}_lca.tsv"
    output:
        "reports/mlca/mlca2tsv/{my_experiment}.tsv"
    params:
        outdir = "reports/mlca/mlca2tsv/{my_experiment}",
        indir = "reports/mlca/",
        rerep = "results/rereplicated/"
    shell:
        "python scripts/mlca-tsv.py -i {params.indir} -r {params.rerep} -o {params.outdir}"



#------------------------------------------------------
# Converting mlca tsv to krona friendly input
#-------------------------------------------------------
rule mlca2kronatext:
    input:
        "results/mlca/{library}/{sample}_lca.tsv"
    output:
        "results/mlca/krona_input/{library}/{sample}_lca.kronatext.tsv"
    shell:
        "sed -r 's/=/\t/g' {input} | tail -n +2 | cut -f 2,5-10 > {output}"


#-----------------------------------------------------
# Krona, interactive html graphics of taxonomic diversity
#-----------------------------------------------------
rule mlca_to_krona:
    conda:
        "../envs/environment.yaml"
    input:
        "results/mlca/krona_input/{library}/{sample}_lca.kronatext.tsv"
    output:
        "reports/krona/mlca/{library}/{sample}.html"
    shell:
        "ktImportText {input} -o {output}"

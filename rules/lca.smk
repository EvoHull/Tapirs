# ==================================================
# LCA: LOWEST COMMON ANCESTOR ANALYSES
# ==================================================

configfile: "config.yaml"

# --------------------------------------------------
# MLCA, majority lowest common ancestor
# --------------------------------------------------

rule mlca:
    conda:
        "../envs/environment.yaml"
    input:
        "results/blasttax/{library}/{sample}.tax.tsv"
    output:
        "results/mlca/{library}/{sample}.lca.tsv"
    params:
        # out = "results/mlca/{library}/{sample}.lca.tsv", # redundant
        bitscore = "10",   # -b blast hit bitscore upper threshold
        identity = "100",  # -id percent identity
        coverage = "60",   # -cov percentage coverage
        majority = "100",  # -m majority percent, 100 is all hits share taxonomy
        hits = "1"  # -hits minimum number of hits, default = 2, 1 is true LCA just takes top hit
    shell:
        "python scripts/mlca.py \
        -i {input} \
        -o {output} \
        -b {params.bitscore} \
        -id {params.identity} \
        -cov {params.coverage} \
        -m {params.majority} \
        -hits {params.hits} \
        "
# GS - The mlca script needs changing because of an error. The final species name is output seperated by a tab and not an undderscore or space; is this something you can fix?

# -------------------------------------------------------
# mlca to tsv
# -------------------------------------------------------

rule mlca_to_tsv:
    conda:
        "../envs/environment.yaml"
    input:
        expand("results/mlca/{library}/{sample}.lca.tsv")
    output:
        "reports/{config[my_experiment]}.tsv"
    params:
        indir = "results/mlca/",
        rerep = "results/rereplicated"  # syntax
    shell:
        "python scripts/mlca-tsv.py -i {params.indir} -r {params.rerep} -o {output}"

# -------------------------------------------------------
# blca, Bayesian lowest common ancestor
# -------------------------------------------------------

rule blca:
    conda:
        "../envs/environment.yaml"
    input:
        query = "results/06_dechimera/{library}/{sample}.nonchimera.fasta",
        database = config["blast_db"]
    output:
        "results/blca/{library}/{sample}.blca.out"
    script:
        "scripts/2.blca_main.py \
        -i {input.query} \
        --db {input.database} \
        -o {output} \
        "
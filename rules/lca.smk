# ==================================================
# LCA: LOWEST COMMON ANCESTOR ANALYSES
# ==================================================

configfile: "config.yaml"

# --------------------------------------------------
# MLCA, majority lowest common ancestor
# --------------------------------------------------

rule mlca:
    input:
        blast = "results/blast_tax/{LIBRARIES}/{SAMPLES}.blast.tax.tsv"
    output:
        lca = "results/mlca/{LIBRARIES}/{SAMPLES}.lca.tsv"
    params:
        bitscore = config['MLCA_bitscore'],
        identity = config['MLCA_identity'],
        coverage = config['MCLA_coverage'],
        majority = config['MLCA_majority'],
        min_hits = config['MLCA_hits']
    script:
        "../scripts/mlca.py"


# MLCA TO TSV

rule mlca_to_tsv:
    input:
        lca = expand("results/mlca/{real_combos}.lca.tsv", real_combos = real_combos),
        rerep = expand("results/09_rereplicated/{real_combos}.rerep.fasta", real_combos = real_combos)
    output:
        tsv = "results/" + config['my_experiment'] + "blast" + config['MLCA_identity'] + ".tsv"
    script:
        "../scripts/mlca_to_tsv.py"

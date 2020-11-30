# ==================================================
# LCA: LOWEST COMMON ANCESTOR ANALYSES
# ==================================================

configfile: "config.yaml"

# --------------------------------------------------
# MLCA, majority lowest common ancestor
# --------------------------------------------------

rule mlca:
    input:
        blast = "results/blast_tax/{SAMPLES}.blast.tax.tsv"
    output:
        lca = "results/mlca/{SAMPLES}.lca.tsv"
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
        lca = expand("results/mlca/{sample}.lca.tsv", sample = SAMPLES),
        rerep = expand("results/09_rereplicated/{sample}.rerep.fasta", sample = SAMPLES)
    output:
        tsv = "results/" + config['my_experiment'] + "blast" + config['MLCA_identity'] + ".tsv"
    script:
        "../scripts/mlca_to_tsv.py"

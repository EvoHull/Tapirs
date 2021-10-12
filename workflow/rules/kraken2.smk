# ==============================================================================
# KRAKEN 2 ANALYSIS
# ==============================================================================

# ------------------------------------------------------------------------------
# KRAKEN 2 TAXONOMIC ASSIGNMENT
# ------------------------------------------------------------------------------

rule kraken2:
    input:
        reads = "results/09_rereplicated/{LIBRARIES}/{SAMPLES}.rerep.fasta"
    output:
        reports = "results/kraken2/reports/{LIBRARIES}/{SAMPLES}.txt",
        outputs = "results/kraken2/outputs/{LIBRARIES}/{SAMPLES}.krk"
    run:
        if len([1 for line in open(input.reads)]) > 0:
            shell("kraken2 --db {config[kraken2_db]} {input.reads} \
            --threads {config[kraken2_threads]} \
            --confidence {config[kraken2_confidence]} \
            --report {output.reports} \
            --output {output.outputs}")
        else:
            open(output.reports, 'w').close()
            open(output.outputs, 'w').close()

# ------------------------------------------------------------------------------
# TAXONOMY TO KRAKEN 2
# ------------------------------------------------------------------------------

rule taxonomy_to_kraken2:
    input:
        config['taxdump'] + '/names.dmp',
        kraken2 = "results/kraken2/outputs/{LIBRARIES}/{SAMPLES}.krk"
    params:
        taxdump = config['taxdump']
    output:
        kraken2_tax = "results/kraken2_tax/{LIBRARIES}/{SAMPLES}.krk.tax.tsv"
    script:
        "../scripts/taxonomy_to_kraken2.py"

# ------------------------------------------------------------------------------
# KRAKEN 2 TO TSV
# ------------------------------------------------------------------------------

rule kraken2_to_tsv:
    input:
        lca = expand("results/kraken2_tax/{real_combos}.krk.tax.tsv", real_combos = real_combos),
        rerep = expand("results/09_rereplicated/{real_combos}.rerep.fasta", real_combos = real_combos)
    params:
        lowest_rank = config['lowest_taxonomic_rank'],
        highest_rank = config['highest_taxonomic_rank']
    output:
        tsv = "results/" + config['my_experiment'] + "_kraken2_conf" + str(config['kraken2_confidence']).split('.')[1] + "_" + config['cluster_method'] + ".tsv"
    script:
        "../scripts/mlca_to_tsv.py"

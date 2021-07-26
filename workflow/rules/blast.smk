# ==================================================
# BLAST ANALYSIS
# ==================================================

# configfile: "config.yaml"

# --------------------------------------------------
# blastn, sequence similarity search
# --------------------------------------------------

rule blastn:
    conda:
        config['conda']
    input:
        query = "results/08_dechimera/{LIBRARIES}/{SAMPLES}.nc.fasta",
    output:
        blast = "results/blast/{LIBRARIES}/{SAMPLES}.blast.tsv"
    params:
        outformat = "'6 qseqid stitle sacc staxids pident qcovs evalue bitscore'",
        evalue = float(config["BLAST_min_evalue"])
    shell:
        "blastn \
        -query {input.query} \
        -db {config[blast_db]} \
        -outfmt {params.outformat} \
        -perc_identity {config[BLAST_min_perc_ident]} \
        -evalue {params.evalue} \
        -max_target_seqs {config[BLAST_max_target_seqs]} \
        -out {output.blast}"

# ------------------------------------------------------------------------------
# TAXONOMY TO BLAST

rule taxonomy_to_blast:
    conda:
        config['conda']
    input:
        taxdump = config['taxdump'],
        blast = "results/blast/{LIBRARIES}/{SAMPLES}.blast.tsv"
    output:
        blast_tax = "results/blast_tax/{LIBRARIES}/{SAMPLES}.blast.tax.tsv"
    script:
        "../scripts/taxonomy_to_blast.py"

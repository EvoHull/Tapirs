# ==================================================
# BLAST ANALYSIS
# ==================================================

configfile: "config.yaml"

# --------------------------------------------------
# blastn, sequence similarity search
# --------------------------------------------------

rule blastn:
    conda:
        "../envs/environment.yaml"
    input:
        query = "results/06_dechimera/{library}/{sample}.nonchimera.fasta"
    params:
        db_dir = config["blast_db"],
        outformat = "'6 qseqid stitle sacc staxids pident qcovs evalue bitscore'"
    output:
        "results/blast/{library}/{sample}.blast.tsv"
    threads:
        6
    shell:
        "blastn \
        -db {params.db_dir} \
        -num_threads {threads} \
        -outfmt {params.outformat} \
        -perc_identity {config[BLAST_min_perc_ident]} \
        -evalue {config[BLAST_min_evalue]} \
        -max_target_seqs {config[BLAST_max_target_seqs]} \
        -query {input.query} \
        -out {output} \
        "

# -----------------------------------------------------
# tax_to_blast, adds taxonomy column to blast output
# -----------------------------------------------------

rule add_taxonomy_to_blast:
    conda:
        "../envs/environment.yaml"
    input:
        blast_out = "results/blast/{library}/{sample}.blast.tsv",
        ranked_lineage = "data/databases/new_taxdump/rankedlineage.dmp"
    output:
        "results/blasttax/{library}/{sample}.tax.tsv",
    params:
    #     taxdump = "data/databases/new_taxdump/rankedlineage.dmp"
    #     indir = "results/blast",
    #     outdir = "results/blasttax"
    shell:
        "python scripts/tax_to_blast.py -i results/blast -o results/blasttax -lin {input.ranked_lineage}"

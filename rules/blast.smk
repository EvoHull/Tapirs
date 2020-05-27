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
        query = expand("results/06_dechimera/{sample.library}/{sample.sample}.nonchimera.fasta",
            sample=sample.reset_index().itertuples()),
    params:
        outformat = "'6 qseqid stitle sacc staxids pident qcovs evalue bitscore'"
    output:
        "results/blast/{library}/{sample}.blast.tsv"
    threads:
        6
    shell:
        "blastn \
        -db {config[blast_db]} \
        -query {input.query} \
        -num_threads {threads} \
        -outfmt {params.outformat} \
        -perc_identity {config[BLAST_min_perc_ident]} \
        -evalue {config[BLAST_min_evalue]} \
        -max_target_seqs {config[BLAST_max_target_seqs]} \
        -out {output} \
        "

# -----------------------------------------------------
# tax_to_blast, adds taxonomy column to blast output
# snakemakification still in progress
# -----------------------------------------------------

rule add_taxonomy_to_blast:
    conda:
        "../envs/environment.yaml"
    input:
        blast_out = expand("results/blast/{sample.library}/{sample.sample}.blast.tsv",
            sample=sample.reset_index().itertuples()),
        ranked_lineage = "data/databases/new_taxdump/rankedlineage.dmp"
    output:
        "results/blasttax/{library}/{sample}.tax.tsv",
    params:
    #     taxdump = "data/databases/new_taxdump/rankedlineage.dmp"
    #     indir = "results/blast",
    #     outdir = "results/blasttax"
    script:
        "scripts/tax_to_blast.py -i results/blast -o results/blasttax -lin {input.ranked_lineage}"

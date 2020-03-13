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
        #db = "nt", #specify in environment.yaml
        query = "results/03_denoised/{library}/{sample}_nc.fasta"
    params:
        db_dir = directory(config["blast_db"]), # database directory
        outformat = "'6 qseqid stitle sacc staxids pident qcovs evalue bitscore'"
    output:
        "results/blast/{library}/{sample}_blast.out"
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
        blast_out = "results/blast/{library}/{sample}_blast.out",
        ranked_lineage = "data/databases/new_taxdump/rankedlineage.dmp"
    output:
        blast_taxonomy = "results/blast/{library}/{sample}_tax.tsv"
    shell:
        "python scripts/tax_to_blast.py -i {input.blast_out} -o {output.blast_taxonomy} -lin {input.ranked_lineage}"

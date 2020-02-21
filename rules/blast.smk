#####   BLASTing   #####
#-----------------------------------------------------
# blastn, sequence similarity search
#-----------------------------------------------------

configfile: "config.yaml"

rule blastn:
    #message: "executing blast analsyis of sequences against database {input.database}"
    conda:
        "../envs/tapirs.yaml"
    input:
        #db = "nt", #specify in environment.yaml
        query = "results/03_denoised/{library}/nc_{sample}.fasta"
    params:
        db_dir = directory("data/databases/12S_full/12s_full"), # database directory
        descriptions = "50", # return maximum of 50 hits
        outformat = "'6 qseqid stitle sacc staxids pident qcovs evalue bitscore'",
        min_perc_ident = "100", # this needs to be 100%
        min_evalue = "1e-20"
    output:
        "results/blast/{library}/{sample}_blast.out"
    threads:
        6
    shell:
        "blastn \
        -db {params.db_dir} \
        -num_threads {threads} \
        -outfmt {params.outformat} \
        -perc_identity {params.min_perc_ident} \
        -evalue {params.min_evalue} \
        -max_target_seqs {params.descriptions} \
        -query {input.query} \
        -out {output} \
        "


#-----------------------------------------------------
# tax_to_blast, adds taxonomy in a column to blast output
#-----------------------------------------------------

rule tax_to_blast:
    conda:
        "../envs/tapirs.yaml"
    input:
        blast_out = "results/blast/{library}/{sample}_blast.out",
        ranked_lineage = "data/databases/new_taxdump/rankedlineage.dmp"
    output:
        blast_taxonomy = "results/blast/{library}/{sample}_tax.tsv"
    shell:
        "python scripts/tax_to_blast.py -i {input.blast_out} -o {output.blast_taxonomy} -lin {input.ranked_lineage}"

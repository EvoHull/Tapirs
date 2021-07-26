# ==================================================
# VSEARCH
# dereplicate, denoise, dechimera, rereplicate
# ==================================================

# configfile: "config.yaml"

# ------------------------------------------------------------------------------
# DEREPLICATIE READS

rule vsearch_dereplicate:
    conda:
        config['conda']
    input:
        fa = "results/05_forward_merged/{LIBRARIES}/{SAMPLES}.fasta"
    output:
        derep = "results/06_dereplicated/{LIBRARIES}/{SAMPLES}.derep.fasta",
        cluster_file = "results/08_clusters/derep/{LIBRARIES}/{SAMPLES}.cluster.tsv"
    shell:
        "vsearch --derep_fulllength {input.fa} \
        --sizeout \
        --minuniquesize {config[VSEARCH_minuniqsize]} \
        --output {output.derep} \
        --uc {output.cluster_file}"

# ------------------------------------------------------------------------------
# CLUSTER READS

rule vsearch_cluster:
    conda:
        config['conda']
    input:
        derep = "results/06_dereplicated/{LIBRARIES}/{SAMPLES}.derep.fasta"
    output:
        cluster = "results/07_clustered/{LIBRARIES}/{SAMPLES}.cluster.fasta",
        cluster_file = "results/08_clusters/cluster/{LIBRARIES}/{SAMPLES}.cluster.tsv"
    shell:
        "vsearch --cluster_fast {input.derep} \
        --sizein --sizeout \
        --query_cov {config[VSEARCH_query_cov]} \
        --id {config[VSEARCH_cluster_id]} \
        --strand both \
        --centroids {output.cluster} \
        --uc {output.cluster_file}"

# ------------------------------------------------------------------------------
# CHIMERA DETECTION

rule vsearch_dechimera:
    conda:
        config['conda']
    input:
        cluster = "results/07_clustered/{LIBRARIES}/{SAMPLES}.cluster.fasta"
    output:
        nonchimeras = "results/08_dechimera/{LIBRARIES}/{SAMPLES}.nc.fasta",
        chimeras = "results/08_dechimera/{LIBRARIES}/{SAMPLES}.chimera.fasta"
    shell:
        "vsearch --uchime_ref {input.cluster} \
        --db {config[dechim_blast_db]} \
        --chimeras {output.chimeras} \
        --borderline {output.chimeras} \
        --mindiffs {config[VSEARCH_mindiffs]} \
        --mindiv {config[VSEARCH_mindiv]} \
        --nonchimeras {output.nonchimeras}"

# ------------------------------------------------------------------------------
# REREPLICATE READS

rule vsearch_rereplicate:
    conda:
        config['conda']
    input:
        nonchimeras = "results/08_dechimera/{LIBRARIES}/{SAMPLES}.nc.fasta"
    output:
        rerep = "results/09_rereplicated/{LIBRARIES}/{SAMPLES}.rerep.fasta",
    shell:
        "vsearch --rereplicate {input.nonchimeras} \
        --output {output.rerep}"

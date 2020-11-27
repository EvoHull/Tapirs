# ==================================================
# VSEARCH
# dereplicate, cluster, dechimera
# ==================================================

configfile: "config.yaml"

# ------------------------------------------------------------------------------
# DEREPLICATIE READS

rule vsearch_dereplicate:
    input:
        fa = "results/05_forward_merged/{SAMPLES}.fasta"
    output:
        derep = "results/06_dereplicated/{SAMPLES}.derep.fasta",
    shell:
        "vsearch --derep_fulllength {input.fa} \
        --sizeout \
        --minuniquesize {config[VSEARCH_minuniqsize]} \
        --output {output.derep} "

# -----------------------------------------------------
# VSEARCH DENOISE

rule vsearch_denoise:
    input:
        "results/06_dereplicated/{SAMPLES}.derep.fasta",
    output:
        seqs = "results/07_denoised/{SAMPLES}.denoise.fasta",
        denoise_results = "reports/vsearch/{SAMPLES}.denoise-report.txt"
    shell:
        "vsearch \
        --cluster_unoise {input} \
        --sizein \
        --sizeout \
        --minsize {config[VSEARCH_minsize]} \
        --unoise_alpha {config[VSEARCH_unoise_alpha]} \
        --id {config[VSEARCH_unoise_id]} \
        --centroids {output.seqs} \
        --uc {output.denoise_results}"

# ------------------------------------------------------------------------------
# CLUSTER READS

# rule vsearch_cluster:
#     input:
#         derep = "results/06_dereplicated/{SAMPLES}.derep.fasta"
#     output:
#         cluster = "results/07_clustered/{SAMPLES}.cluster.fasta",
#         cluster_file = "results/08_clusters/{SAMPLES}.cluster.txt"
#     shell:
#         "vsearch --cluster_fast {input.derep} \
#         --sizein --sizeout \
#         --query_cov {config[VSEARCH_query_cov]} \
#         --id {config[VSEARCH_cluster_id]} \
#         --strand both \
#         --centroids {output.cluster} \
#         --uc {output.cluster_file}"

# ------------------------------------------------------------------------------
# CHIMERA DETECTION

rule vsearch_dechimera:
    input:
        cluster = "results/07_denoised/{SAMPLES}.denoise.fasta"
    output:
        nonchimeras = "results/08_dechimera/{SAMPLES}.nc.fasta",
        chimeras = "results/08_dechimera/{SAMPLES}.chimera.fasta"
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
    input:
        nonchimeras = "results/08_dechimera/{SAMPLES}.nc.fasta"
    output:
        rerep = "results/09_rereplicated/{SAMPLES}.rerep.fasta",
    shell:
        "vsearch --rereplicate {input.nonchimeras} \
        --output {output.rerep}"

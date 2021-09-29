# ==================================================
# VSEARCH
# dereplicate, cluster/denoise, chimera detection, rereplicate
# ==================================================

# ------------------------------------------------------------------------------
# DEREPLICATIE READS
# ------------------------------------------------------------------------------

rule vsearch_dereplicate:
    conda:
        config['conda']
    input:
        fa = "results/05_fasta/{LIBRARIES}/{SAMPLES}.fasta"
    output:
        derep = "results/06_dereplicated/{LIBRARIES}/{SAMPLES}.derep.fasta",
        cluster_file = "results/08_clusters/derep/{LIBRARIES}/{SAMPLES}.cluster.tsv"
    shell:
        "vsearch \
        --derep_fulllength {input.fa} \
        --sizeout \
        --minuniquesize {config[VSEARCH_minuniqsize]} \
        --output {output.derep} \
        --fasta_width 0 \
        --uc {output.cluster_file}"

# ------------------------------------------------------------------------------
# CLUSTER READS
# ------------------------------------------------------------------------------

rule vsearch_cluster:
    conda:
        config['conda']
    input:
        derep = "results/06_dereplicated/{LIBRARIES}/{SAMPLES}.derep.fasta"  if config['cluster_method'] == "cluster" else []
    output:
        cluster = "results/07_clustered/{LIBRARIES}/{SAMPLES}.cluster.fasta",
        cluster_file = "results/08_clusters/cluster/{LIBRARIES}/{SAMPLES}.cluster.tsv"
    shell:
        "vsearch \
        --cluster_fast {input.derep} \
        --sizein --sizeout \
        --query_cov {config[VSEARCH_query_cov]} \
        --id {config[VSEARCH_cluster_id]} \
        --strand both \
        --centroids {output.cluster} \
        --fasta_width 0 \
        --uc {output.cluster_file}"

# ------------------------------------------------------------------------------
# DENOISE READS
# ------------------------------------------------------------------------------

rule vsearch_denoise:
    conda:
        config['conda']
    input:
        "results/06_dereplicated/{LIBRARIES}/{SAMPLES}.derep.fasta" if config['cluster_method'] == "denoise" else []
    output:
        seqs = "results/07_denoised/{LIBRARIES}/{SAMPLES}.denoise.fasta",
        denoise_results = "results/08_clusters/denoise/{LIBRARIES}/{SAMPLES}.denoise.tsv"
    shell:
        "vsearch \
        --cluster_unoise {input} \
        --sizein --sizeout \
        --minsize {config[VSEARCH_minsize]} \
        --unoise_alpha {config[VSEARCH_unoise_alpha]} \
        --id {config[VSEARCH_unoise_id]} \
        --centroids {output.seqs} \
        --fasta_width 0 \
        --uc {output.denoise_results}"

# ------------------------------------------------------------------------------
# CHIMERA DETECTION
# ------------------------------------------------------------------------------

if config['chimera_detection'] == "ref":

    rule vsearch_uchime_ref:
        conda:
            config['conda']
        input:
            cluster = "results/07_clustered/{LIBRARIES}/{SAMPLES}.cluster.fasta" if config['cluster_method'] == "cluster" else "results/07_denoised/{LIBRARIES}/{SAMPLES}.denoise.fasta"
        output:
            nonchimeras = "results/08_dechimera/{LIBRARIES}/{SAMPLES}.nc.fasta",
            chimeras = "results/08_dechimera/{LIBRARIES}/{SAMPLES}.chimera.fasta"
        shell:
            "vsearch \
            --uchime_ref {input.cluster} \
            --db {config[dechim_blast_db]} \
            --chimeras {output.chimeras} \
            --borderline {output.chimeras} \
            --mindiffs {config[VSEARCH_mindiffs]} \
            --mindiv {config[VSEARCH_mindiv]} \
            --fasta_width 0 \
            --nonchimeras {output.nonchimeras}"

elif config['chimera_detection'] == "denovo":

    rule vsearch_uchime3_denovo:
        conda:
            config['conda']
        input:
            cluster = "results/07_clustered/{LIBRARIES}/{SAMPLES}.cluster.fasta" if config['cluster_method'] == "cluster" else "results/07_denoised/{LIBRARIES}/{SAMPLES}.denoise.fasta"
        output:
            nonchimeras = "results/08_dechimera/{LIBRARIES}/{SAMPLES}.nc.fasta",
            chimeras = "results/08_dechimera/{LIBRARIES}/{SAMPLES}.chimera.fasta"
        shell:
            "vsearch \
            --uchime3_denovo {input.cluster} \
            --abskew {config[VSEARCH_abskew]} \
            --chimeras {output.chimeras} \
            --borderline {output.chimeras} \
            --fasta_width 0 \
            --nonchimeras {output.nonchimeras}"

# ------------------------------------------------------------------------------
# REREPLICATE READS
# ------------------------------------------------------------------------------

rule vsearch_rereplicate:
    conda:
        config['conda']
    input:
        nonchimeras = "results/08_dechimera/{LIBRARIES}/{SAMPLES}.nc.fasta"
    output:
        rerep = "results/09_rereplicated/{LIBRARIES}/{SAMPLES}.rerep.fasta",
    shell:
        "vsearch \
        --rereplicate {input.nonchimeras} \
        --fasta_width 0 \
        --output {output.rerep}"

# ==================================================
# ESTABLISH TAXONOMY DATA DIR
# ==================================================

configfile: "config.yaml"

# --------------------------------------------------
# Get taxonomy dump data from NCBI
# --------------------------------------------------

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

rule get_new_taxdump:
    conda:
        "../envs/environment.yaml"
    threads:
        8
    input:
        FTP.remote("ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip")
    output:
        "data/databases/new_taxdump/rankedlineage.dmp"
    params:
        "data/databases/new_taxdump"
    shell:
        "unzip {input} -d {params}"

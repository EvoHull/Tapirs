# ==================================================
# ESTABLISH TAXONOMY DATA DIR
# ==================================================

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
configfile: "config/config.yaml"

# --------------------------------------------------
# Get taxonomy dump data from NCBI
# --------------------------------------------------

FTP = FTPRemoteProvider()

rule get_new_taxdump:
    conda:
        "../envs/environment.yaml"
    threads:
        8
    input:
        FTP.remote("ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip")
    output:
        "resources/databases/new_taxdump/rankedlineage.dmp"
    params:
        "resources/databases/new_taxdump"
    shell:
        "unzip {input} -d {params}"

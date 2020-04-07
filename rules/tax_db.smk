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
    input:
        FTP.remote("ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip")
    output:
        directory("data/databases/new_taxdump")
    shell:
        "unzip {input} -d {output}"

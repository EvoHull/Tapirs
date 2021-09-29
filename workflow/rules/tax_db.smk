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
        config['conda']
    input:
        FTP.remote("ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip")
    output:
        config['taxdump'] + "/names.dmp"
    params:
        config['taxdump']
    shell:
        "unzip {input} -d {params}"

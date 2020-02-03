#-----------------------------------------------------
# metaclese-blast
# ---------
# A metabarcoding workflow to make taxon assignments
# using blast
#-----------------------------------------------------

# do I need to specify config.yaml here? Or does it inherit from top level snakemake file?
configfile: "config.yaml" # specify a configuration file

# get the sequence files into a list
SAMPLES, = glob_wildcards("results/seqkit/{sample}")


#-----------------------------------------------------
# target rule
#-----------------------------------------------------
rule all_blast:
    input:
        expand("results/seqkit/{sample}", sample=SAMPLES),
        expand("results/blast/config[expt_name].out", sample=SAMPLES),
        expand("results/basta/{sample}.basta_LCA.out", sample=SAMPLES),
        expand("results/basta/{sample}.basta_LCA.out.biom", sample=SAMPLES)


#-----------------------------------------------------
# make blast database
# this is done once and not strictly part of the workflow
#-----------------------------------------------------
# rule blast_db:
#     input:
#         "data/database/database.fas.gz"
#     output:
#         destination=directory("results/blast/database"),
#         name="my_db_name"
#     shell: # is this fastest way to process gz data?
#         "gunzip -c {input} | makeblastdb -in - -dbtype nucl -title my_db -out {output.name}"


#-----------------------------------------------------
# blastn - get top hit
#-----------------------------------------------------
rule blast:
    input:
        database="config[blast_db]", # blast db is in data/blast as specified in config file
        query="results/seqkit/{sample}.extendedFrags.fas"
    params:
        descriptions="50", # return maximum of 50 hits
        outformat="'6 qseqid stitle pident evalue'"
    output:
        "results/blast/{config[expt_name]}.out" # experiment name from a config.yaml file
    shell:
        "blastn -db {input.database} -outfmt {params.outformat} -max_target_seqs {params.descriptions} -query {input.query} -out {output}"

#-----------------------------------------------------
# BASTA - Last Comomon Ancestor analysis of blast
#-----------------------------------------------------
rule basta:
    input:
        "results/blast/{sample}_blastn.out" #fix this
        # file of blast tabular -outfmt 6 from above
    params: # need to check --help to see how to flag these at shell
        percentID="100%", # must be shared by all, 90% = taxonomy shared by 90% hits, try this
        minhits="3", # must have at least 3 hits, else ignored
        evalue= "10e‐20",#  try 1e‐8
        length="100"
    conda:
        "envs/basta_env.yaml" # Basta needs a python 2.7 environment
    output:
        "results/basta/{sample}.basta_LCA.out"
    shell:
        "./bin/basta sequence {input} {output} gb -p {params.percentID} -e {params.evalue} -l {params.length}"
#        "./bin/basta multiple INPUT_DIRECTORY OUTPUT_FILE MAPPING_FILE_TYPE" from manual

#-----------------------------------------------------
# BASTA to BIOM - BASTA output tsv needs to go into BIOM
# uses BIOM-convert
#-----------------------------------------------------
rule basta_BIOM:
    input:
        "results/basta/basta_LCA.out"
    params:
        json="json",
        hdf5="hdf5"
    output:
        "results/basta/{sample}.basta_LCA.out.biom"
    shell:
        "biom convert -i {input} -o {output} --table-type='OTU table' --to-{params.json}"

# biom convert -i table.txt -o table.from_txt_json.biom --table-type="OTU table" --to-json
# biom convert -i table.txt -o table.from_txt_hdf5.biom --table-type="OTU table" --to-hdf5

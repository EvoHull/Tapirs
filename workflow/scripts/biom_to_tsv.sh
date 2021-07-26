# This is the command to convert biom files into tsv files - biom is already loaded in conda in this environment

# -i denotes the input file and -o denotes the output file
biom convert -i testing_inverts_tapirs.kraken2.biom -o table.from_biom.tsv --to-tsv

# Need to run this after to include sample names in the tsv table
sed -i "/OTU/s/^#//g" table.from_biom.tsv

#Then rscript can be run for analyses - will also need to load R in the environment
R vegan_analyses_test.R

import pandas as pd
outfile = 'test.csv'
# library = pd.read_csv("libraries.csv", sep='\t', header=0, index_col="library")
sample = pd.read_csv("samples.tsv", sep='\t', header=0, index_col="library", dtype=str)
sample.to_csv(outfile)
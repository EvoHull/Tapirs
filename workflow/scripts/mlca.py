# Majority Lowest Common Ancestor (mlca)
# --------------------------------------
# caluclation of mlca from blast output
# part of Tapirs metabarcoding workflow

# # import libraries
import pandas as pd
import numpy as np

# input data in blast format
infile = snakemake.input.blast
# results file, tabular format
outfile = snakemake.output.lca

identity = float(snakemake.params.identity)
majority = float(snakemake.params.majority)
coverage = float(snakemake.params.coverage)
min_hits = float(snakemake.params.min_hits)
bit_threshold = float(snakemake.params.bitscore)
prop = 1 - (bit_threshold / 100)

unid_tax = ('\t'.join([str(x) for x in ['unidentified']*7]))
taxonomy = ('kingdom','phylum','class','order','family','genus','species')
header = 'query\tlca_rank\tlca_otu\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tmethod\n'

with open(infile, 'r') as file:
    count=0
    for line in file.readlines():
        count += 1
    if count == 0:
        with open(outfile, 'w') as lca_out:
            lca_out.write(header)
    else:
        df = pd.read_csv(infile, sep = '\t', header=None)
        df = df[df[4] >= identity]
        df = df[df[5] >= coverage]
        seq_ids = sorted(set(df[0]))

        with open(outfile, 'w') as lca_out:
            lca_out.write(header)
            for seq_id in seq_ids:
                lca = []
                seq = df[df[0] == seq_id]
                seq = seq[seq[7] >= max(seq[7]) * prop]
                taxon = seq[8].drop_duplicates()
                taxon = taxon.str.split("/", expand = True)

                for col in taxon.columns:
                    taxa = np.unique(taxon[col], return_counts = True)
                    if max(taxa[1]) / sum(taxa[1]) >= majority / 100:
                        lca.append(taxa[0][taxa[1] == max(taxa[1])])

                lca_tax = []
                if len(lca) == 0:
                    lca.append('unidentified')
                for i in lca:
                    lca_tax.append(i[0])
                if len(lca_tax) < 7:
                    lca_tax = lca_tax + (['unidentified'] * (7 - len(lca_tax)))
                else:
                    lca_tax = lca_tax
                lca_tax = ('\t'.join([str(x) for x in lca_tax]))

                otu_id = str(lca[-1][0])

                if len(seq[col]) >= min_hits and len(seq[col]) > 1:
                    lca_out.write('%s\t%s\t%s\t%s\tlca\n' %(seq_id, taxonomy[len(lca)-1], otu_id, lca_tax))
                elif len(seq[col]) == 1 and min_hits == 1:
                    lca_out.write('%s\t%s\t%s\t%s\tsingle_hit\n' %(seq_id, taxonomy[len(lca)-1], otu_id, lca_tax))
                else:
                    lca_out.write('%s\tunidentified\tunidentified\t%s\tunidentified\n' %(seq_id, unid_tax))

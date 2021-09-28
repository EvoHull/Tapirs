# ------------------------------------------------------------------------------
# MAJORITY LOWEST COMMON ANCESTOR (mlca)
# caluclation of mlca from blast output
# ------------------------------------------------------------------------------

import pandas as pd
import numpy as np
from itertools import dropwhile

infile = snakemake.input.blast
outfile = snakemake.output.lca

identity = float(snakemake.params.identity)
majority = float(snakemake.params.majority)
coverage = float(snakemake.params.coverage)
min_hits = float(snakemake.params.min_hits)
bit_threshold = float(snakemake.params.bitscore)
prop = 1 - (bit_threshold / 100)

unid_tax = ('\t'.join([str(x) for x in ['unidentified']*7]))
taxonomy = ('domain','phylum','class','order','family','genus','species')
header = 'query\ttax_rank\totu_id\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\tmethod\n'

# remove trailing unknowns from lca taxonomy list
def strip_unknown(l):
        l = list(dropwhile(lambda x: x == 'unknown', l[::-1]))
        return l[::-1]

with open(infile, 'r') as file:
    count = len([1 for line in open(infile)])  # check if file is not an empty input file
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

                if len(seq) >= min_hits:
                    taxon = seq[8].drop_duplicates()
                    taxon = taxon.str.split("/", expand = True)
                    for col in taxon.columns:
                        taxa = np.unique(taxon[col], return_counts = True)
                        if max(taxa[1]) / sum(taxa[1]) >= majority / 100:
                            lca.append(taxa[0][taxa[1] == max(taxa[1])])

                    if len(seq) >= min_hits and len(seq) > 1:
                        lca_method = 'lca'
                    elif len(seq) == 1 and min_hits == 1:
                        lca_method = 'single_hit'

                    lca_tax = []
                    lca = strip_unknown(lca)
                    if len(lca) > 0:
                        for i in lca:
                            lca_tax.append(i[0])
                        if len(lca_tax) < 7:
                            lca_tax = lca_tax + (['unidentified'] * (7 - len(lca_tax)))
                        lca_tax = ('\t'.join([str(x) for x in lca_tax]))

                        otu_id = str(lca[-1][0])

                        lca_out.write('%s\t%s\t%s\t%s\t%s\n' %(seq_id, taxonomy[len(lca)-1], otu_id, lca_tax, lca_method))

                    else:
                        lca_out.write('%s\tunidentified\tunidentified\t%s\t%s\n' %(seq_id, unid_tax, lca_method))

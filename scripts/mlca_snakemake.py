# Majority Lowest Common Ancestor (mlca)
# --------------------------------------
# caluclation of mlca from blast output
# part of Tapirs metabarcoding workflow

# import libraries
import argparse
import pandas as pd
import numpy as np

# argument parser for variables
# parser = argparse.ArgumentParser(description='')
# parser.add_argument('-i', '--infile', metavar='blast output', dest='infile', type=str,
#             help='input data in blast format', default='', required=True)
# parser.add_argument('-o', '--outfile',metavar='output file', dest='outfile', type=str,
#             help='results file in tabular', required=True)
# parser.add_argument('-b', '--bitscore', metavar='bitscore percentage threshold', dest='bit_threshold', type=str,
#             help='top precentage threshold for hit bitscore', required=True)
# parser.add_argument('-id', metavar='identity', dest='identity', type=str,
#             help='identity threshold', required=True)
# parser.add_argument('-cov', metavar='coverage', dest='coverage', type=str,
#             help='coverage threshold', required=True)
# parser.add_argument('-hits', metavar='minimum hits', dest='min_hits', type=str,
#             help='minimum hits required for a query to be processed, set to 1 will do stuff', default='2', required=False)
# parser.add_argument('-m', metavar='majority', dest='majority', type=str,
#             help='majority percentage hits to be identical at taxonomic rank for assignment', required=True)
# args = parser.parse_args()

# input data in blast format
infile=snakemake.input[0]
# results file, tabular format
outfile=snakemake.output[0]

# definitions
bit_threshold=float(snakemake.params.bitscore)
identity=float(snakemake.params.identity)
coverage=float(snakemake.params.coverage)
majority=float(snakemake.params.majority)
min_hits=float(snakemake.params.hits)

# clarify minimum number of hits
if min_hits=='':
    min_hits=2
else:
    min_hits=min_hits # ???
prop=1-(bit_threshold/100)

# unidentified taxa
unid_tax=('\t'.join([str(x) for x in ['unidentified']*7]))
taxonomy=('kingdom','phylum','class','order','family','genus','species')

# write output document
with open(infile,'r') as file:
    count=0
    for line in file.readlines():
        count+=1
    if count == 0:
        with open(outfile,'w') as lca_out:
            lca_out.write('query\tlca_rank\tlca_otu\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tmethod\n')
    else:
        df=pd.read_csv(infile,sep='\t', header=None)
        df=df[df[4]>=identity]
        df=df[df[5]>=coverage]
        seq_ids=set(df[0])

        with open(outfile,'w') as lca_out:
            lca_out.write('query\tlca_rank\tlca_otu\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tmethod\n')
            for seq_id in seq_ids:
                lca=[]
                seq=df[df[0] == seq_id]
        # already defined prop above?
                prop=1-(bit_threshold/100)
                seq=seq[seq[7]>=max(seq[7])*prop]
                taxon=seq[8].str.split("/", expand = True)
                for col in taxon.columns:
                    taxa=np.unique(taxon[col],return_counts=True)
                    if max(taxa[1])/len(taxon[col]) >= majority/100:
                        lca.append(taxa[0][taxa[1]==max(taxa[1])])
                # denote as 'unidentified'
                lca_tax=[]
                if len(lca)==0:
                    lca.append('unidentified')
                for i in lca:
                    lca_tax.append(i[0])
                if len(lca_tax)<7:
                    lca_tax=lca_tax+(['unidentified']*(7-len(lca_tax)))
                else:
                    lca_tax=lca_tax
                lca_tax=('\t'.join([str(x) for x in lca_tax]))

                if len(lca)==7:
                    otu_id=('_'.join([str(lca[5][0]),str(lca[6][0])]))
                elif len(lca)==6:
                    otu_id=('_'.join([str(lca[5][0]),'spp.']))
                else:
                    otu_id=str(lca[-1][0])

                if len(seq[col])>=min_hits and len(seq[col]) > 1:
                    lca_out.write('%s\t%s\t%s\t%s\tlca\n' %(seq_id,taxonomy[len(lca)-1],otu_id,lca_tax))
                elif len(seq[col]) == 1 and min_hits == 1:
                    lca_out.write('%s\t%s\t%s\t%s\tsingle_hit\n' %(seq_id,taxonomy[len(lca)-1],otu_id,lca_tax))
                else:
                    lca_out.write('%s\tunidentified\tunidentified\t%s\tunidentified\n' %(seq_id,unid_tax))

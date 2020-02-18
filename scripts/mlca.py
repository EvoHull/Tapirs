import argparse
import pandas as pd
import numpy as np

# usage:
# python tax_to_blast.py -i blast_out/BLE04_blast.tsv -o blast_out/BLE04_tax.tsv -lin new_taxdump/rankedlineage.dmp
#       tax_to_blast.py adds taxonomy in a column to blast output
#       rankedlineage.dmp is genbank taxonomy

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--infile', metavar='blast output', dest='infile', type=str,
            help='input data in blast format', default='', required=True)
parser.add_argument('-o', '--outfile',metavar='output file', dest='outfile', type=str,
            help='results file in tabular', required=True)
parser.add_argument('-b', '--bitscore', metavar='bitscore percentage threshold', dest='bit_threshold', type=str,
            help='top precentage threshold for hit bitscore', required=True)
parser.add_argument('-id', metavar='identity', dest='identity', type=str,
            help='identity treshold', required=True)
parser.add_argument('-cov', metavar='coverage', dest='coverage', type=str,
            help='coverage treshold', required=True)
parser.add_argument('-hits', metavar='minimum hits', dest='min_hits', type=str,
            help='minimum hits required for a query to be processed, set to 1 will do stuff', default='2', required=False)
parser.add_argument('-m', metavar='majority', dest='majority', type=str,
            help='majority percentage hits to be identical at taxonomic rank for assignment', required=True)
args = parser.parse_args()

infile=args.infile
outfile=args.outfile

identity=float(args.identity)
majority=float(args.majority)
coverage=float(args.coverage)
if float(args.min_hits)=='':
    min_hits=2
else:
    min_hits=float(args.min_hits)
bit_threshold=float(args.bit_threshold)

df=pd.read_csv(args.infile,sep='\t', header=None)
df[df[4]>=float(args.identity)]
df[df[4]>=float(args.coverage)]
seq_ids=set(df[0])
prop=1-(float(args.bit_threshold)/100)

taxonomy=('kingdom','phylum','class','order','family','genus','species')

df=pd.read_csv(infile,sep='\t', header=None)
df[df[4]>=identity]
df[df[4]>=coverage]
seq_ids=set(df[0])

with open(outfile,'w') as lca_out:
    lca_out.write('query\tlca_rank\tlca_taxon\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tmethod\n')
    for seq_id in seq_ids:
        lca=[]
        seq=df[df[0] == seq_id]
        prop=1-(bit_threshold/100)
        seq=seq[seq[7]>=max(seq[7])*prop]
        taxon=seq[8].str.split("/", expand = True)
        for col in taxon.columns:
            taxa=np.unique(taxon[col],return_counts=True)
            if max(taxa[1])/len(taxon[col]) >= majority/100:
                lca.append(taxa[0][taxa[1]==max(taxa[1])])

        lca_tax=[]
        for i in lca:
            lca_tax.append(i[0])
        if len(lca_tax)<7:
            lca_tax=lca_tax+(['unidentified']*(7-len(lca_tax)))
        else:
            lca_tax=lca_tax
        lca_tax=('\t'.join([str(x) for x in lca_tax]))
        unid_tax=('\t'.join([str(x) for x in ['unidentified']*7]))
        if len(seq[col])>=min_hits and len(seq[col]) > 1:
            lca_out.write('%s\t%s\t%s\t%s\tlca\n' %(seq_id,taxonomy[len(lca)-1],str(lca[-1][0]),lca_tax))
        elif len(seq[col]) == 1 and min_hits == 1:
            lca_out.write('%s\t%s\t%s\t%s\tsingle_hit\n' %(seq_id,taxonomy[len(lca)-1],str(lca[-1][0]),lca_tax))
        else:
            lca_out.write('%s\tunidentified\tunidentified\t%s\tunidentified\n' %(seq_id,unid_tax))

import argparse
import glob
import os

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--indir', metavar='blast input directory', dest='indir', type=str,
            help='directory containing libraries with blast outputs', default='', required=True)
parser.add_argument('-o', '--outdir',metavar='output directory', dest='outdir', type=str,
            help='directory for taxonomy assigned libraries', required=True)
parser.add_argument('-lin', '--lineage',metavar='ranked lineage', dest='ranked', type=str,
            help='results file in tsv', required=True)
args = parser.parse_args()

def reference_taxonomy():
    taxonomyDict = {}
    with open(args.ranked,'r') as rankedlineage:
        for tax in rankedlineage:
            tax = tax.split('|')
            taxonid = tax[0]
            species = tax[1].strip().replace(' ','_') if tax[1].strip() else 'unknown_species'
            genus = tax[3].strip() if tax[3].strip() else 'unknown'
            family = tax[4].strip() if tax[4].strip() else 'unknown'
            order = tax[5].strip() if tax[5].strip() else 'unknown'
            classe = tax[6].strip() if tax[6].strip() else 'unknown'
            phylum = tax[7].strip() if tax[7].strip() else 'unknown'
            kingdom = tax[8].strip() if tax[8].strip() else 'unknown'
            superkingdom = tax[9].strip() if tax[9].strip() else 'unknown'
            taxonomyDict[str(tax[0].strip())] = {'species':species, 'genus':genus, 'family':family, 'order':order, 'class':classe, 'phylum':phylum, 'kingdom':kingdom,'superkingdom':superkingdom}
    return taxonomyDict

taxonomyDict=reference_taxonomy()

libraries=sorted(glob.glob(args.indir+'/*'))
for library in libraries:
    library_dir=args.outdir+'/'+library.split('/')[-1]
    if not os.path.exists(library_dir):
        os.mkdir(library_dir)
    files=sorted(glob.glob(library+'/*_blast.tsv'))
    
    for file in files:
        sample=file.split('/')[-1].split('_')[0]
        with open(file,'r') as blasthits, open(library_dir+'/'+sample+'_tax.tsv', 'w') as output:
            for hit in blasthits:
                taxid = hit.split('\t')[3]
                if taxid == 'N/A':
                    output.write(hit.strip()+'\tunknown/unknown/unknown/unknown/unknown/unknown/unknown\n')
                else:
                    output.write(hit.strip()+'\t'+taxonomyDict[taxid]['superkingdom'] \
                                 +'/'+taxonomyDict[taxid]['phylum'] \
                                 +'/'+taxonomyDict[taxid]['class'] \
                                 +'/'+taxonomyDict[taxid]['order'] \
                                 +'/'+taxonomyDict[taxid]['family'] \
                                 +'/'+taxonomyDict[taxid]['genus'] \
                                 +'/'+taxonomyDict[taxid]['species'].split('_')[1]+'\n')

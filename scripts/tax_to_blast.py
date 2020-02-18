import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--infile', metavar='blast output', dest='infile', type=str,
            help='input data in blast format', required=True)
parser.add_argument('-o', '--outfile',metavar='output file', dest='outfile', type=str,
            help='results file in tsv', required=True)
parser.add_argument('-lin', '--lineage',metavar='ranked lineage', dest='ranked', type=str,
            help='results file in tsv', required=True)
args = parser.parse_args()

taxonomyDict = {}
with open(args.ranked,'r') as rankedlineage:
    for tax in rankedlineage:
        tax = tax.split('|')
        taxonid = tax[0]
        species = tax[1].strip().replace(' ','_') if tax[1].strip() else 'unknown_species'
        genus = tax[3].strip() if tax[3].strip() else 'unknown_genus'
        family = tax[4].strip() if tax[4].strip() else 'unknown_family'
        order = tax[5].strip() if tax[5].strip() else 'unknown_order'
        classe = tax[6].strip() if tax[6].strip() else 'unknown_class'
        phylum = tax[7].strip() if tax[7].strip() else 'unknown_phylum'
        kingdom = tax[8].strip() if tax[8].strip() else 'unknown_kingdom'
        superkingdom = tax[9].strip() if tax[9].strip() else 'unknown_superkingdom'
        taxonomyDict[str(tax[0].strip())] = {'species':species, 'genus':genus, 'family':family, 'order':order, 'class':classe, 'phylum':phylum, 'kingdom':kingdom,'superkingdom':superkingdom}


with open(args.infile,'r') as blasthits, open(args.outfile, 'w') as output:
    for hit in blasthits:
        taxid = hit.split('\t')[3]
        if taxid == 'N/A':
            output.write(hit.strip()+'\tunknown_kingdom/unknown_phylum/unknown_class/unknown_order/unknown_family/unknown_genus/unknown_species\n')
        else:
            output.write(hit.strip()+'\t'+taxonomyDict[taxid]['superkingdom']+'/'+taxonomyDict[taxid]['phylum']+'/'+taxonomyDict[taxid]['class']+'/'+taxonomyDict[taxid]['order']+'/'+taxonomyDict[taxid]['family']+'/'+taxonomyDict[taxid]['genus']+'/'+taxonomyDict[taxid]['species']+'\n')

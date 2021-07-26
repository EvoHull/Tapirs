
# version 2
# this uses merged.dmp to correct taxids from rankedlineage.dmp
# to use this, direct snakemake to the 'new_taxdump' directory rather than rankedlineage.dmp
# do this in the config file at TAXDUMP:

def reference_taxonomy():
    taxonomyDict = {}
    with open(snakemake.input.taxdump + '/rankedlineage.dmp', 'r') as rankedlineage:
        for tax in rankedlineage:
            tax = tax.split('|')
            taxonid = tax[0]
            species = tax[1].strip().replace(
                ' ', '_') if tax[1].strip() else 'unknown'
            genus = tax[3].strip() if tax[3].strip() else 'unknown'
            family = tax[4].strip() if tax[4].strip() else 'unknown'
            order = tax[5].strip() if tax[5].strip() else 'unknown'
            classe = tax[6].strip() if tax[6].strip() else 'unknown'
            phylum = tax[7].strip() if tax[7].strip() else 'unknown'
            kingdom = tax[8].strip() if tax[8].strip() else 'unknown'
            superkingdom = tax[9].strip() if tax[9].strip() else 'unknown'
            taxonomyDict[str(tax[0].strip())] = {'species': species, 'genus': genus, 'family': family,
                                                 'order': order, 'class': classe, 'phylum': phylum, 'kingdom': kingdom, 'superkingdom': superkingdom}
    return taxonomyDict

def merged_taxonomy():
    mergedDict = {}
    with open(snakemake.input.taxdump + '/merged.dmp', 'r') as merged:
        for taxid in merged:
            taxid = taxid.split('|')
            mergedDict[taxid[0].strip()] = taxid[1].strip()
    return mergedDict

taxonomyDict = reference_taxonomy()
mergedDict = merged_taxonomy()

with open(snakemake.input.blast, 'r') as blasthits, open(snakemake.output.blast_tax, 'w') as output:
    for hit in blasthits:
        taxid = hit.split('\t')[3].split(';')[0]
        if taxid == 'N/A':
            output.write(hit.strip() + '\tunknown/unknown/unknown/unknown/unknown/unknown/unknown\n')
        else:
            try:
                superkingdom = taxonomyDict[taxid]['superkingdom']
            except KeyError:
                taxid = mergedDict[taxid]
                superkingdom = taxonomyDict[taxid]['superkingdom']

            if superkingdom != 'unknown':
                output.write(hit.strip() + '\t' + taxonomyDict[taxid]['superkingdom']
                             + '/' + taxonomyDict[taxid]['phylum']
                             + '/' + taxonomyDict[taxid]['class']
                             + '/' + taxonomyDict[taxid]['order']
                             + '/' + taxonomyDict[taxid]['family']
                             + '/' + taxonomyDict[taxid]['genus']
                             + '/' + '_'.join(taxonomyDict[taxid]['species'].split('_')[:2]) + '\n')
            else:
                output.write(hit.strip() + '\t' + taxonomyDict[taxid]['kingdom']
                             + '/' + taxonomyDict[taxid]['phylum']
                             + '/' + taxonomyDict[taxid]['class']
                             + '/' + taxonomyDict[taxid]['order']
                             + '/' + taxonomyDict[taxid]['family']
                             + '/' + taxonomyDict[taxid]['genus']
                             + '/' + '_'.join(taxonomyDict[taxid]['species'].split('_')[:2]) + '\n')

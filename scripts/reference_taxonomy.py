
def reference_taxonomy():
    taxonomyDict = {}
    with open(snakemake.input.taxdump, 'r') as rankedlineage:
        for tax in rankedlineage:
            tax = tax.split('|')
            taxonid = tax[0]
            species = tax[1].strip().replace(
                ' ', '_') if tax[1].strip() else 'unknown_species'
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


reference_taxonomy()
#snakemake.output.taxonomyDict = taxonomyDict

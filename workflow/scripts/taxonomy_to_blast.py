# adapted from taxopy core.py functions (https://github.com/apcamargo/taxopy)
# includes merged taxid dictionary for taxid corrections adapted from Simple-LCA (https://github.com/sdwfrost/Simple-LCA)

nodes_dmp = snakemake.params.taxdump + '/nodes.dmp'
names_dmp = snakemake.params.taxdump + '/names.dmp'
merged_dmp = snakemake.params.taxdump + '/merged.dmp'

blast = snakemake.input.blast
blast_tax = snakemake.output.blast_tax

# adapted from taxopy core.py
def import_nodes():
    taxid2parent = {}
    taxid2rank = {}
    with open(nodes_dmp, 'r') as file:
        for line in file:
            line = line.split('\t')
            taxid = int(line[0])
            parent = int(line[2])
            rank = line[4]
            taxid2parent[taxid] = parent
            taxid2rank[taxid] = rank
    return taxid2parent, taxid2rank

# adapted from taxopy core.py
def import_names():
    taxid2name = {}
    with open(names_dmp, 'r') as file:
        for line in file:
            line = line.split('\t')
            if line[6] == 'scientific name':
                taxid = int(line[0])
                name = line[2]
                taxid2name[taxid] = name
    return taxid2name

# adapted from Simple-LCA
def merged_taxonomy():
    merged_dict = {}
    with open(merged_dmp, 'r') as merged:
        for taxid in merged:
            taxid = taxid.split('|')
            merged_dict[taxid[0].strip()] = taxid[1].strip()
    return merged_dict

# adapted from taxopy core.py
def get_lineage(taxid):
    lineage = []
    current_taxid = taxid
    lineage.append(current_taxid)
    while taxid2parent[current_taxid] != current_taxid:
        current_taxid = taxid2parent[current_taxid]
        lineage.append(current_taxid)
    return lineage

# adapted from taxopy core.py
def rank_name_dictionary(taxid):
    lineage = get_lineage(taxid)
    rank_name_dictionary = {}
    for taxid in lineage:
        rank = taxid2rank[taxid]
        if rank != 'no rank':
            rank_name_dictionary[rank] = taxid2name[taxid]
    return rank_name_dictionary

def taxonomy_string(taxid):
    # taxonomic ranks
    ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    # generating taxonomic rank dictionary from input taxid
    taxdict = rank_name_dictionary(taxid)
    # dropping genus_sp. and species without genus
    if 'species' in taxdict:
        taxdict['species'] = taxdict['species'].split('/')[0]
        taxdict['species'] = '_'.join(taxdict['species'].split(' ')[0:2])
        if '_sp.' in taxdict['species']:
            del(taxdict['species'])
    if 'species' in taxdict and 'genus' not in taxdict:
        del(taxdict['species'])
    # fill in missing taxonomic ranks with 'unknown'
    tax_ranks = []
    for rank in ranks:
        if rank in list(taxdict.keys()):
            tax_ranks.append(taxdict[rank])
        else:
            tax_ranks.append('unknown')
    return tax_ranks

# if input file is not empty generate taxonomy dictionaries
if len([1 for line in open(blast)]) > 0:
    taxid2parent, taxid2rank = import_nodes()
    taxid2name = import_names()
    merged_dict = merged_taxonomy()

    # write taxonomy string to output file using taxonomy rank names from taxid
    with open(blast, 'r') as blasthits, open(blast_tax, 'w') as output:
        for hit in blasthits:
            taxid = hit.split('\t')[3]
            if taxid == 'N/A':
                output.write(hit.strip() + '\tunknown/unknown/unknown/unknown/unknown/unknown/unknown\n')
            else:
                if ';' not in str(taxid):
                    try:
                        output.write(hit.strip() + '\t' + '/'.join(taxonomy_string(int(taxid))) + '\n')
                    except KeyError:
                        taxid = merged_dict[taxid]
                        output.write(hit.strip() + '\t' + '/'.join(taxonomy_string(int(taxid))) + '\n')

# if input file is empty write empty file to output
else:
    open(blast_tax, 'w').close()

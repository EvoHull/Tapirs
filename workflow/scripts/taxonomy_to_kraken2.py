# adapted from taxopy core.py functions (https://github.com/apcamargo/taxopy)
# includes merged taxid dictionary for taxid corrections adapted from Simple-LCA (https://github.com/sdwfrost/Simple-LCA)

from itertools import dropwhile

nodes_dmp = snakemake.params.taxdump + '/nodes.dmp'
names_dmp = snakemake.params.taxdump + '/names.dmp'
merged_dmp = snakemake.params.taxdump + '/merged.dmp'

kraken = snakemake.input.kraken2
kraken_tax = snakemake.output.kraken2_tax

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

def strip_unknown(l):
        l = list(dropwhile(lambda x: x == 'unknown', l[::-1]))
        return l[::-1]

# taxonomy level names and output file header
taxonomy = ('domain','phylum','class','order','family','genus','species')  # taxonomy level names used
header = 'query\ttax_rank\totu_id\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\n'  # header for output file

# if input file is not empty generate taxonomy dictionaries
if len([1 for line in open(kraken)]) > 0:
    taxid2parent, taxid2rank = import_nodes()
    taxid2name = import_names()
    merged_dict = merged_taxonomy()

    # write taxonomy string to output file using taxonomy rank names from taxid
    with open(kraken, 'r') as kraken_in, open(kraken_tax, 'w') as output:
        output.write(header)
        for hit in kraken_in:
            hit = hit.split('\t')
            if hit[0] == 'C':  # select only classified
                taxid = hit[2]
                if taxid != '1':  # ignore 'root' level assignments
                    tax_str = strip_unknown(taxonomy_string(int(taxid)))  # trim trailing unknowns on taxonomy string
                    tax_rank = taxonomy[len(tax_str)-1]  # get lowest assigned taxonomic rank
                    otu_id = tax_str[-1]  # get lowest assigned taxonomic rank name
                    tax_str = tax_str + (['unidentified'] * (7 - len(tax_str)))  # fill taxonomy to length with unidentified
                    tax_str = ('\t'.join([str(x) for x in tax_str]))
                    output.write('%s\t%s\t%s\t%s\n' %(hit[1].strip(), tax_rank, otu_id, tax_str))

# if input file is empty write empty file to output
else:
    with open(kraken_tax, 'w') as output:
        output.write(header)

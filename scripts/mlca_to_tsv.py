# ----------------------------------------------------
# tsv output from mlca
# ----------------------------------------------------

# import argparse
import pandas as pd
import glob

# argparse can be replaced by snakemake
# ----------------------------------------------------
# parser = argparse.ArgumentParser(description='')
# parser.add_argument('-i', '--indir', metavar='mlca output', dest='indir', type=str,
#             help='directory containing mlca outputs', default='', required=True)
# parser.add_argument('-r', '--rerep',metavar='rerep output', dest='rerep', type=str,
#             help='directory containing vsearch rereplicated outputs', required=True)
# parser.add_argument('-o', '--outfile',metavar='output file', dest='outfile', type=str,
#             help='output tsv file', required=True)
# args = parser.parse_args()


final_out = pd.DataFrame()
taxonomy = pd.DataFrame()
pd.set_option('max_colwidth', 400)

outfile = snakemake.output.tsv

samples = snakemake.input.lca
samples.sort()
for sample in samples:
    name = sample.split('.')[0]
    sample_name = name.split('/')[-1]
    name = [name.split('/')[-2], name.split('/')[-1]]
    name = '/'.join(name)

    with open(sample, 'r') as file:  # open sample file
        count = 0
        for line in file.readlines():  # count number of lines
            count += 1
            if count > 1:
                df = pd.read_csv(file.name, sep='\t',
                                 header=None, skiprows=1)
                df[0] = pd.to_numeric(
                    df[0].str.split(r";size=", expand=True)[1])
            else:
                df = pd.DataFrame([([0])+(['unidentified']*10)])

    species_list = sorted(set(df[2][df[2] != 'unidentified']))
    tsv_index = species_list
    tsv_index.append('unassigned')
    index = pd.Index(tsv_index)

    dfob = pd.DataFrame(columns=[sample_name], index=index)

    df = df.replace('unidentified', '')
    tdf = pd.concat([df[df.columns[2:10]]], axis=1, sort=False).drop_duplicates(
        2).sort_values(2).set_index(2)

    for species in species_list:
        reads = df[0][df[2] == species].sum()
        dfob.loc[species] = reads

    assigned_reads = dfob[sample_name].sum()

    rerep = 'results/09_rereplicated/' + name + '.rerep.fasta'

    total_reads = len([1 for line in open(rerep)
                       if line.startswith('>')])

    unassigned_reads = total_reads-assigned_reads
    dfob.loc['unassigned'] = unassigned_reads

    taxonomy = (pd.concat([taxonomy, tdf], axis=0,
                          sort=True)).drop_duplicates().sort_index()

    final_out = (pd.concat([final_out, dfob], axis=1, sort=False)).fillna(
        0).astype(int).sort_index()

tfob = pd.DataFrame(columns=['taxonomy'], index=final_out.index)
for otu in final_out.index:
    if otu != 'unassigned':
        ranks = taxonomy[taxonomy.index == otu].to_string(
            index=False, header=False).split('  ')
        tax_add = ('; '.join(['k__'+ranks[0].replace(' ', ''),
                              'p__'+ranks[1],
                              'c__'+ranks[2],
                              'o__'+ranks[3],
                              'f__'+ranks[4],
                              'g__'+ranks[5],
                              's__'+ranks[6]]))
        tfob.loc[otu] = tax_add
    else:
        tfob.loc[otu] = 'u__unassigned'

final_out = (pd.concat([final_out, tfob], axis=1, sort=False)).fillna(0).sort_index()
final_out.index.name = '#OTU_ID'
final_out = final_out.drop('unassigned', axis=0).append(final_out.loc[['unassigned'], :])
final_out.to_csv(outfile, sep='\t', index=True, header=True)

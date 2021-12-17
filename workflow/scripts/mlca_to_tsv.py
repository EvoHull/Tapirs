# ------------------------------------------------------------------------------
# tsv output from lca and tax_to_kraken2 outputs
# part of Tapirs metabarcoding workflow

# ------------------------------------------------------------------------------

# import libraries
import pandas as pd
import glob

final_out = pd.DataFrame()  # create empty final dataframe to populate
taxonomy = pd.DataFrame()  # create empty taxonomy dataframe to populate
pd.set_option('max_colwidth', 400)  # allow longer strings in dataframes

outfile = snakemake.output.tsv
samples = snakemake.input.lca
minimum_rank = snakemake.params.lowest_rank
maximum_rank = snakemake.params.highest_rank

samples = sorted(samples)  # sort all input files alphabetically

for sample in samples:
    name = sample.split('.')[0]  # drop extensions from file path
    sample_name = name.split('/')[-1]  # sample name from file path
    name = [name.split('/')[-2], name.split('/')[-1]]  # library and sample names from file path
    name = '/'.join(name)  # make library/sample as string

    with open(sample, 'r') as file:  # open sample file
        count = len([1 for line in open(sample)])  # check if file is not an empty input file
        if count > 1:
            df = pd.read_csv(file.name, sep = '\t', header = None, skiprows = 1)  # dataframe from mlca file and drop headers
            if ';size=' in str(df[0]):  # using vsearch cluster size extension as read count
                df[0] = pd.to_numeric(df[0].str.split(r";size=", expand = True)[1])  # cut out read name and keep cluster size
            else:  # for kraken2 output
                df[0] = 1  # each read is single
        else:  #if file is an empty input file
            df = pd.DataFrame([([0]) + (['unidentified'] * 10)])  # creat dataframe with '0' unidentified read mlca output

    ranks = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    min_rank = int(ranks.index(minimum_rank))
    max_rank = int(ranks.index(maximum_rank))
    ranks_upper = ranks[max_rank:]  # list maximum rank required and all ranks below
    ranks_lower = ranks[min_rank + 1:]  # list all ranks below minimum rank required

    for i in range(min_rank + 4, 10):  # change all ranks below minimum rank required to unidentified
        df[i] = 'unidentified'

    df = df[df[1].isin(ranks_upper)]  # trim dataframe of assignments higher than maximum rank required

    rank_list = list(df[1])  # list existing assigned ranks
    for n, i in enumerate(rank_list):  # loop through assigned ranks and...
        if i in ranks_lower:
            rank_list[n] = ranks[min_rank]  # ...change assigned ranks lower to minimum required rank
    df[1] = rank_list

    lis = list(df[min_rank + 3])  # list assigned otu minimum required rank names
    lit = list(df[2])  # list otu assigned rank names
    for i in range(len(lis)):
        if lis[i] != 'unidentified':
            lit[i] = lis[i]  # change assigned otu names to minimum required otu rank names if lower than required minimum rank
    df[2] = lit

    species_list = sorted(set(df[2][df[2] != 'unidentified']))  # create otu list
    tsv_index = species_list
    tsv_index.append('unassigned')  # add unassigned to otu list
    index = pd.Index(tsv_index)  # and make otu list a dataframe index

    dfob = pd.DataFrame(columns=[sample_name], index=index)

    df = df.replace('unidentified', '')
    tdf = pd.concat([df[df.columns[2:10]]], axis = 1, sort = False).drop_duplicates(2).sort_values(2).set_index(2)

    for species in species_list:  # get read counts per otu
        reads = df[0][df[2] == species].sum()
        dfob.loc[species] = reads  # place reads per otu in dataframe

    assigned_reads = dfob[sample_name].sum()  # get total assigned reads

    rerep = 'results/09_rereplicated/' + name + '.rerep.fasta'  # grab rereplicated fasta for total reads

    total_reads = len([1 for line in open(str(rerep)) if line.startswith('>')])  # total reads calculation

    unassigned_reads = total_reads - assigned_reads  # unassigned reads calculation
    dfob.loc['unassigned'] = unassigned_reads  # add unassigned reads to dataframe

    taxonomy = (pd.concat([taxonomy, tdf], axis = 0, sort = True)).drop_duplicates().sort_index()

    final_out = (pd.concat([final_out, dfob], axis = 1, sort = False)).fillna(0).astype(int).sort_index()

final_out = final_out.reindex(sorted(final_out.columns), axis=1)  # order columns alphabetically

tfob = pd.DataFrame(columns = ['taxonomy'], index = final_out.index)  # create, sort and write taxonomy strings
for otu in final_out.index:
    if otu != 'unassigned':
        ranks = taxonomy[taxonomy.index == otu].to_string(index = False, header = False).split(' ')
        tax_add = ('; '.join(['d__' + ranks[0].replace(' ', ''),
                              'p__' + ranks[1],
                              'c__' + ranks[2],
                              'o__' + ranks[3],
                              'f__' + ranks[4],
                              'g__' + ranks[5],
                              's__' + ranks[6]]))
        tax_add_ranks = [y for y in [x.split('__')[1] for x in tax_add.split('; ')] if y]
        tfob.loc[otu] = '; '.join(tax_add.split('; ')[0:len(tax_add_ranks)])
    else:
        tfob.loc[otu] = 'u__unassigned'

final_out = (pd.concat([final_out, tfob], axis = 1, sort = False)).fillna(0).sort_index()  # add taxonomy column to final dataframe

# order final dataframe rows by taxonomy alphabetically and by ranks length: species -> domain
str_len = [len(y) for y in [x.split(';') for x in final_out['taxonomy']]]  # lengths of taxonomy strings
final_out['str_len'] = str_len  # adding column of taxonomy string length
final_out = final_out.sort_values(['str_len','taxonomy'], ascending = [False, True])  # sort by taxonomy alphabetically and by ranks length
final_out.drop('str_len', inplace = True, axis = 1)  # drop taxonomy string length column

final_out.index.name = '#OTU_ID'
final_out = final_out.drop('unassigned', axis = 0).append(final_out.loc[['unassigned'], :])  # unassigned reads as last row
final_out.to_csv(outfile, sep = '\t', index = True, header = True)  # write the output tsv file

import argparse
import pandas as pd
import numpy as np
import glob

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--indir', metavar='mlca output', dest='indir', type=str,
            help='directory containing mlca outputs to be combined', default='', required=True)
parser.add_argument('-o', '--outfile',metavar='output file', dest='outfile', type=str,
            help='output tsv file', required=True)
args = parser.parse_args()

final_out=pd.DataFrame()
taxonomy=pd.DataFrame()
pd.set_option('max_colwidth',400)

files=sorted(glob.glob(args.indir+'*_lca.tsv'))
for file in files:
    file_name=file.split('/')[-1].split('_')[0]
    
    with open(file,'r') as file:
        count=0
        for line in file.readlines():
            count+=1
            if count > 1:
                df=pd.read_csv(file.name, sep='\t', header=None,skiprows=1)
                df[0]=pd.to_numeric(df[0].str.split(r"size=", expand=True)[1])
            else:
                df=pd.DataFrame([([0])+(['unidentified']*10)])

    species_list=sorted(set(df[2][df[2]!= 'unidentified']))
    tsv_index=species_list
    tsv_index.append('unassigned')
    index=pd.Index(tsv_index)
    
    dfob=pd.DataFrame(columns=[file_name],index=index)
    
    df=df.replace('unidentified','')
    tdf=pd.concat([df[df.columns[2:10]]], axis=1, sort=False).drop_duplicates(2).sort_values(2).set_index(2)
    
    for species in species_list:
        reads=df[0][df[2]==species].sum()
        dfob.loc[species]=reads
    
    assigned_reads=dfob[file_name].sum()
    
    total_reads=len([1 for line in open('6_denoise_uc/rerep/EA01/'+file_name+'_rerep.fasta') \
                     if line.startswith('>')])
    
    unassigned_reads=total_reads-assigned_reads
    dfob.loc['unassigned']=unassigned_reads
    
    taxonomy=(pd.concat([taxonomy, tdf], axis=0, sort=True)).drop_duplicates().sort_index()
    
    final_out=(pd.concat([final_out, dfob], axis=1, sort=False)).fillna(0).astype(int).sort_index()

tfob=pd.DataFrame(columns=['taxonomy'],index=final_out.index)
for otu in final_out.index:
    if otu != 'unassigned':
        ranks=taxonomy[taxonomy.index==otu].to_string(index=False,header=False).split('  ')
        tax_add=('; '.join(['k__'+ranks[0].replace(' ',''), \
                            'p__'+ranks[1], \
                            'c__'+ranks[2], \
                            'o__'+ranks[3], \
                            'f__'+ranks[4], \
                            'g__'+ranks[5], \
                            's__'+ranks[6]]))
        tfob.loc[otu]=tax_add
    else:
        tfob.loc[otu]='u__unassigned'
        
final_out=(pd.concat([final_out, tfob], axis=1, sort=False)).fillna(0).sort_index()
final_out.index.name='#OTU_ID'
final_out=final_out.drop('unassigned', axis=0).append(final_out.loc[['unassigned'],:])
final_out.to_csv(args.outfile+'.tsv', sep='\t',index=True, header=True)

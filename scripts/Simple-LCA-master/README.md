# Simple-LCA
A tool to simply find the lowest common ancestor from blast results with taxonids. This approach is partly based on MEGAN's (Huson et al., 2007) LCA method.
## Getting Started
### Prerequisites
* Python 2.7
### Installing
Clone this repository
```
git clone https://github.com/naturalis/Simple-LCA
```
### Download NCBI taxonomy files
This script makes use of the new_taxdump files from the NCBI ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/
```
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
```
## Usage
### Quick start
The script uses a specific BLAST output, below an example of the BLAST command
```
blastn -query [input] -out [output] -db nt -max_target_seqs 100 -outfmt '6 qseqid stitle sacc staxids pident qcovs evalue bitscore'
```
Before finding the lca you need to add taxonomy to the blast output
```
python add_taxonomy.py -i [blast file] -t rankedlineage.dmp -m merged.dmp -o [output]
```
Now the taxonomy is added, it is time the find the lowest common ancestor
```
python lca.py -i [blast with taxonomy] -o [output] -b 8 -id 80 -cov 80 -t yes -tid 99 -tcov 100 -fh "environmental" -flh "unknown" 
```
### Parameters
...

## Source
Huson, D. H., Auch, A. F., Qi, J., & Schuster, S. C. (2007). MEGAN analysis of metagenomic data. Genome Research, 17(3), 377–386. http://doi.org/10.1101/gr.5969107

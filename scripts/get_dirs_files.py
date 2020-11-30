# script to write directory and sample names to samples.tsv file
# -------------------------------------------------------------
# entries should correspond to library and sample names
# if sequence files do not end in 'fastq.gz' or 'fq.gz'
# the endings should be ammended below

import os
import pandas as pd

# location of the data directories, must end in /
startdir = "data/01_demultiplexed/"

with open('samples_with_duplicates.tsv', 'w') as outfile:  # create output file samples.tsv
    # outfile.write('library' + '\t' + 'sample' + '\n')  # add column headers
    # iterate across files from startdir down
    for root, dirs, files in os.walk(startdir):
        for file in files:
            if file.endswith(('fastq.gz', 'fq.gz')):  # specify file endings here
                path = os.path.join(root, file)
                head_tail = os.path.split(path)  # split into path & file
                # remove startdir from path leaving directory name 'dirs'
                dirs = head_tail[0].replace(startdir, '')
                filetail = head_tail[1]  # get tail, ie filename
                # get sample name from filename losing '.R1.fastq.gz'
                # split filename on . as list, get first item
                sample_name = filetail.rsplit('.', -1)[0]
                # make list of directories and sample names
                dflist = [dirs, sample_name]
                dfstr = '\t'.join(dflist)  # convert to string with tabs
                # write to output file with newline after dir tab sample
                outfile.write(dfstr + '\n')

# specify input and output files
dupfile = "samples_with_duplicates.tsv"
samples_output = "samples.tsv"

# read in the tsv file containing duplicate names from F and R files
df = pd.read_csv(dupfile, sep="\t")
# remove duplicate rows
df.drop_duplicates(subset=None, inplace=True)
# add column names
# df.columns = ['library', 'sample']
# Write the results to a new tsv, no (index) row names
df.to_csv(samples_output, sep="\t", index=False, header=False)

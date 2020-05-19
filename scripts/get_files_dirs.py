# script to write directory and sample names to samples.tsv file
# -------------------------------------------------------------
# entries should correspond to library and sample names
# if sequence files do not end in 'fastq.gz' or 'fq.gz
# they should be ammended below

import os

startdir = "data/01_demultiplexed/"  # location of the data directories

with open('samples.tsv', 'w') as outfile:  # create output file samples.tsv
    outfile.write('library' + '\t' + 'sample' + '\n')  # add column headers
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
                sample_name = filetail.rsplit('.', -1)[0]  # split filename on . as list, get first item
                dflist = [dirs, sample_name]  # make list of directories and sample names
                dfstr = '\t'.join(dflist)  # convert to string with tabs
                outfile.write(dfstr + '\n')  # write to output file with newline after dir tab sample

# script to print dir and file names to tsv file
import os

startdir = "data/01_demultiplexed/"
with open('samples.tsv', 'w') as outfile:
    outfile.write('Library' + '\t' + 'Sample' + '\n') # add column headers
    for root, dirs, files in os.walk(startdir):  # iterate across files from startdir down
        for file in files:
            if file.endswith(('fastq.gz', 'fq.gz')): # specify file endings
                path= os.path.join(root, file)
                head_tail = os.path.split(path)  # split into path & file
                # remove startdir from path leaving directory name 'dirs'
                dirs = head_tail[0].replace(startdir, '')
                filetail = head_tail[1] # get tail, ie file
                dflist = [dirs, filetail] # make list of directories and files
                dfstr = '\t'.join(dflist) # convert to string with tabs
                outfile.write(dfstr + '\n')  # write newline after dir tab file

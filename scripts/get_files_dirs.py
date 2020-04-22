# script to print file and dir names to tsv file
import os

startdir = "data/01_demultiplexed/"
with open('file_tab_py.tsv', 'w') as outfile:
    for root, dirs, files in os.walk(startdir):
        for file in files: # iterate across files from startdir down
            if file.endswith(('fastq.gz', 'fq.gz')): # specify file endings
                path= os.path.join(root, file)
                head_tail = os.path.split(path)  # split into path & file
                # remove startdir from path leaving directory name 'dirs'
                dirs = head_tail[0].replace(startdir, '')
                dflist = [dirs, file] # make list for csv.writer
                sep = '\t'  # use tab as seperator
                dfstr = sep.join(dflist) # convert to sting with tabs
                outfile.write(dfstr + '\n')  # write newline after dir tab file

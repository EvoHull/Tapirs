# save all directory names to a file

import os

basepath = "data/01_demultiplexed/"  # starting location

with os.scandir(basepath) as entries:
    with open('libraries.tsv', 'w') as outfile:  # specify file to write
        for entry in entries:
            if entry.is_dir():  # find directories
                dirnames = [entry.name]  # list
                s1='\n'.join(dirnames)  # string
                outfile.write(s1+'\n')  # write to file with newlines

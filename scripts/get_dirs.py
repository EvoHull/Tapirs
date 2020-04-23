# save all directory names to a file

import os

basepath = "data/01_demultiplexed/"  # starting location

with os.scandir(basepath) as entries:
    with open('libraries.tsv', 'w') as outfile:  # specify file to write
        outfile.write('Library' + '\n') # add Library column header
        for entry in entries:
            if entry.is_dir():  # find directories
                dirnames = [entry.name]  # list
                # write directory names, as strings, to file
                outfile.write(str("".join(dirnames)) + '\n')

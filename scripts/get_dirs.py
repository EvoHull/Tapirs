# save all directory names to libraries.tsv
# -----------------------------------------
# these dirs correspond to sequencing library names

import os

# location of the data directories, must end in /
basepath = "data/01_demultiplexed/"

with os.scandir(basepath) as entries:
    with open('libraries.tsv', 'w') as outfile:  # specify file to write
        outfile.write('library' + '\n')  # add Library column header
        for entry in entries:
            if entry.is_dir():  # find directories
                dirnames = [entry.name]  # list
                # write directory names, as strings, to file
                outfile.write(str("".join(dirnames)) + '\n')

# insert check for duplicate library names

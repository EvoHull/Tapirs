# script to write directory and sample names to samples.tsv file
# -------------------------------------------------------------
# Entries should correspond to sequence file names ending in
# '.fasta' or '.fas' and their containing directory names.
# The file path is also written.

import os

# location of the data directories, must end in /
STARTDIR = "resources/sequences/"

with open("config/samples.tsv", "w") as outfile:  # create output file samples.tsv
    outfile.write('sample_ID' + '\t' + 'File' + '\t' + 'dir_name' + '\t' + 'Path' + '\n')  # add column headers
    # iterate across files from STARTDIR down
    for root, dirs, files in os.walk(STARTDIR):
        for file in files:
            if file.endswith((".fasta", ".fas", ".fa", ".fna")):  # specify file endings here
                path = os.path.join(root, file)
                head_tail = os.path.split(path)  # split into path & file
                # remove STARTDIR from path leaving directory name 'dirs'
                dirs = head_tail[0].replace(STARTDIR, "")
                filetail = head_tail[1]  # get tail, ie filename
                # split filename on . as list, get first item (part before dot=sample_name)
                sample_name = filetail.rsplit(".", -1)[0]
                # make list of sample names, directories, and paths
                dflist = [sample_name, file, dirs, path]
                dfstr = "\t".join(dflist)  # convert to string with tabs
                # write to output file with newline after dir tab sample
                outfile.write(dfstr + "\n")

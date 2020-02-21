import os

with open("files.txt", "w") as a:
    for path, subdirs, files in os.walk(r'../data/01_demultiplexed'):
       for filename in files:
         f = os.path.join(path, filename)
         a.write(str(f) + os.linesep)

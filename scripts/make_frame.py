import os, glob
import pandas as pd

path = r'data/01_demultiplexed/'

all_files = glob.glob(os.path.join(path, "*/*.fastq.gz"))

all_data = pd.DataFrame([])

for path in all_files:
    file = path.replace(".R1.fastq.gz", "")
    file2 = file.replace(".R2.fastq.gz", "")
    entry = file2.replace("data/01_demultiplexed/", "")
    all_data = all_data.append(pd.DataFrame({'library':entry}, index=[0]), ignore_index=True)

    all_data[['Library','Sample']] = all_data.library.str.split("/",expand=True)

all_data = all_data.drop_duplicates()
del all_data["library"]

all_data.to_csv('output.tsv', header = False, index = False, sep = "\t")

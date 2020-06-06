import pandas as pd
configfile: "config.yaml"

# --------------------------------------------------
# Wildcarding library and sample
# --------------------------------------------------

# library = pd.read_table(config["libraries"], index_col="library")
library = pd.read_csv(config["libraries"], sep='\t',
                      header=0, index_col="library")

# sample = pd.read_table(config["samples"], dtype=str)
# sample = pd.read_table(config["samples"], index_col=["library", "sample"], dtype=str)
sample = pd.read_csv(config["samples"], sep='\t',
                     header=0, index_col="library", dtype=str),
# index_col=["library", "sample"],

print("------")
# sample.head()

print("------")

# for value in sample["sample"]:
#     print(value)

print("------")

sample.to_csv("test.csv")





# print(expand("{sample.library}",
#     sample=sample.reset_index().itertuples())),

# sample.index = sample.index.set_levels([i.astype(str) for i in sample.index.levels])

# print(sample.sample) # + "sample" + "\n" + sample.library + "libraries")
# print(library)
# sample['sample']
# library = pd.read_csv(config["libraries"], sep='\t',
#                      header=0, index_col="library")
# # library = pd.read_table(config["libraries"], index_col="library")
# sample = pd.read_csv(config["samples"], sep='\t',
#                      header=0, index_col=["library", "sample"], dtype=str)
# # sample = pd.read_table(config["samples"], index_col=[
# #                        "library", "sample"], dtype=str)
# sample.MultiIndex = sample.MultiIndex.set_levels(
#     [sample.i.astype(str) for i in sample.MultiIndex.levels])

# sample.index = sample.index.set_levels(
#     [i.astype(str) for i in sample.index.levels])
# df.index = df.index.set_levels(df.index.levels[2].astype(int), level=2)
#  so `i.astype(str)` should be `sample.i.astype(str)`

import os  #1
basepath = ‘../data/01_demultiplexed’  #2
with os.scandir(basepath) as entries:
    for entry in entries:
        if entry.is_dir():  #3.1 4.1
            print(entry.name), #3.2ish
            for f_name in os.listdir(‘basepath.entry’):
                if f_name.endswith(‘.fastq’):
                    print(entry.name, ‘/t’ f_name)

#
#
# import os
# with os.scandir('my_directory/') as entries:
#     for entry in entries:
#         print(entry.name)
#
# # -----------------------
# import os
#
# # List all files in a directory using scandir()
# basepath = 'my_directory/'
# with os.scandir(basepath) as entries:
#     for entry in entries:
#         if entry.is_file():
#             print(entry.name)
#
#   # List all subdirectories using scandir()
#   basepath = 'my_directory/'
#   with os.scandir(basepath) as entries:
#       for entry in entries:
#           if entry.is_dir():
#               print(entry.name)
#
#
#   for f_name in os.listdir('some_directory'):
#   ...     if f_name.endswith('.txt'):
#   ...         print(f_name)

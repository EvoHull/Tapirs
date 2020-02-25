#!/bin/bash

echo -e "library\tsample" > samples.tsv
<<<<<<< HEAD
find data/01_demultiplexed/* -name *fastq* | sed -r 's/\//\t/g' |  cut -f 3- | sed -r 's/\./\t/g' | sed -r 's/_/\t/g' | cut -f 1-2 | sort | uniq >> samples.tsv
=======
find data/01_demultiplexed/* -name *fastq* | sed -r 's/\//\t/g' |  cut -f 3- | sed -r 's/\./\t/g' | cut -f 1-2 | sort | uniq >> samples.tsv
>>>>>>> ad581f057e1660cfbb485711b2161dd0273f171c
cut -f 1 samples.tsv | uniq >> libraries.tsv

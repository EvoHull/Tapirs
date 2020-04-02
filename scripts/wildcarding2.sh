#!/bin/bash

echo -e "library\tsample" > samples.tsv
find data/01_demultiplexed/* -name *fastq* | sed -E --regexp-extended 's/\//\t/g' |  cut -f 3- | sed -E --regexp-extended 's/\./\t/g' | sed -E --regexp-extended 's/_/\t/g' | cut -f 1-2 | sort | uniq >> samples.tsv
cut -f 1 samples.tsv | uniq >> libraries.tsv

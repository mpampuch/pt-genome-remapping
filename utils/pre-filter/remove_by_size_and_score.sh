#!/bin/bash

inputPafFile=$1
size=$2
score=$3

# The first positional argument is the input paf file
# The second positional argument is the minimum size of the alignment
# The third positional argument is the minimum score of the alignment

# Example usage:
# ./remove_by_size_and_score.sh /home/kevin/Downloads/chr21.paf 1000 30

awk -v size=$size -v score=$score '{if ($12 >= score && $9 - $8 >= size) print $0}' $inputPafFile 
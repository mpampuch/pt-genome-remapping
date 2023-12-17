#!/bin/bash

featureList=$1

# Create a conda environment with bioawk and activate it
conda activate ${CONDA_ENV}

arr=()
while IFS= read -r line; do
   arr+=("$line")
done <$featureList

ls -d -1 */ > temp_dir_file.txt

dirs="./temp_dir_file.txt"

dir_arr=()
while IFS= read -r dline
do
  dir_arr+=("$dline")
done < "$dirs"

# explaination on how to index arrays https://stackoverflow.com/questions/6723426/looping-over-arrays-printing-both-index-and-value
for i in "${!arr[@]}"
do
    inFileName="../n_check/${dir_arr[$i]}passed_n_check/${arr[$i]}_features_with_no_Ns.fa"
    outFileNameOneFA="./${arr[$i]}_features_with_gt_or_eq_to_20_bp.fa"
    outFileNameOne="./${arr[$i]}_features_with_gt_or_eq_to_20_bp_w_length.txt"
    outFileNameTwoFA="./${arr[$i]}_features_with_lt_20_bp.fa"
    outFileNameTwo="./${arr[$i]}_features_with_lt_20_bp_w_length.txt"

    bioawk -c fastx '{ if(length($seq) >= 20) { print ">"$name; print $seq }}' < $inFileName > $outFileNameOneFA
    bioawk -c fastx '{ print $name, length($seq) }' < $outFileNameOneFA > $outFileNameOne

    bioawk -c fastx '{ if(length($seq) < 20) { print ">"$name; print $seq }}' < $inFileName > $outFileNameTwoFA
    bioawk -c fastx '{ print $name, length($seq) }' < $outFileNameTwoFA > $outFileNameTwo
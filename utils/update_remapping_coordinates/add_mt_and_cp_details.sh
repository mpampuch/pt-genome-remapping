#!/bin/bash

# This script updates the successfully remapped features with their new coordinates
# It takes the following arguments:
# 1. The path to the gffFile merged with the fully filtered BLAST file sorted by chromosome name and then group by chromosome and sorted by start coordinate
inputSortedGFFandBLASTFile=$1
# 2. The path to the gff file containing all the features from the mitochondrial genome (needs to have cleaned up metadata)
inputMtGFFFile=$2
# 3. The path to the gff file containing all the features from the chloroplast genome (needs to have cleaned up metadata)
inputCpGFFFile=$3

# Step 4: Append the mitochondrial genome information
awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' "${inputMtGFFFile}" >>"${inputSortedGFFandBLASTFile}"

# # Step 5: Append the chloroplast genome information
awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' "${inputCpGFFFile}" >>"${inputSortedGFFandBLASTFile}"

# # Step 6: Remove Empty Lines produced by appending
awk 'NF' "${inputSortedGFFandBLASTFile}" >genome_data/new/final_remapped_data_updated_coordinates_sorted_plus_cp_and_mt_minus_original_metadata.gff

# # Step 7: Replace the metadata in the new gff file with the original metadata (perform this in Python because it's much easier)

# # Comment this out if you want to inspect the intermediate files
# rm genome_data/new/*TMP*

#!/bin/bash

# This script updates the successfully remapped features with their new coordinates
# It takes the following arguments:
# 1. The path to the gffFile merged with the fully filtered BLAST file
inputMergedGFFandBLASTFile=$1
# TODO move the follinwg files to difference script
# 2. The path to the gff file containing all the features from the mitochondrial genome (needs to have cleaned up metadata)
# inputMtGFFFile=$2
# 3. The path to the gff file containing all the features from the chloroplast genome (needs to have cleaned up metadata)
# inputCpGFFFile=$3

# Step 1: Update the chromosome name and the coordinates of the features in the new gff file
## Column 1 is the original chromosome name and need to be replaced with the new chromosome name (column 11)
## Column 4 is the original start coordinate and need to be replaced with the new start coordinate (column 18) if the strand is plus and with the new end coordinate (column 19) if the strand is minus
## Column 5 is the original end coordinate and need to be replaced with the new end coordinate (column 19) if the strand is plus and with the new start coordinate (column 18) if the strand is minus
## Column 7 is the original strand and need to be replaced with the new strand (column 22)
awk -F "\t" 'BEGIN {OFS="\t"} {$1=$11; if ($22=="plus") $4=$18; if ($22=="plus") $5=$19; if ($22=="plus") $7="+"; if ($22=="minus") $4=$19; if ($22=="minus") $5=$18; if ($22=="minus") $7="-"; print $0}' "${inputMergedGFFandBLASTFile}" >genome_data/new/final_remapped_data_updated_coordinates_plus_all_blast_data_minus_cp_and_mt_minus_original_metadata.gff

# Step 2: Print only the GFF columns
awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' genome_data/new/final_remapped_data_updated_coordinates_plus_all_blast_data_minus_cp_and_mt_minus_original_metadata.gff >genome_data/new/final_remapped_data_updated_coordinates_minus_cp_and_mt_minus_original_metadata.gff

# Step 3: Sort the new gff file by chromosome name and then group by chromosome and sort by start coordinate (NEED R to do group_by)
# sort -n -t "_" -k 2 ${inputMergedGFFandBLASTFile}_TMP1.gff > ${inputMergedGFFandBLASTFile}_TMP2.gff

# TODO MOVE THE FOLLOWING CODE INTO A NEW FILE
# Step 4: Append the mitochondrial genome information
# awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' ${inputMtGFFFile} >> ${inputMergedGFFandBLASTFile}_TMP2.gff

# # Step 5: Append the chloroplast genome information
# awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' ${inputCpGFFFile} >> ${inputMergedGFFandBLASTFile}_TMP2.gff

# # Step 6: Remove Empty Lines produced by appending
# awk 'NF' ${inputMergedGFFandBLASTFile}_TMP2.gff > genome_data/new/final_remapped_data_MINUS_ORIGINAL_METADATA.gff

# # Step 7: Replace the metadata in the new gff file with the original metadata (perform this in Python because it's much easier)

# # Comment this out if you want to inspect the intermediate files
# rm genome_data/new/*TMP*

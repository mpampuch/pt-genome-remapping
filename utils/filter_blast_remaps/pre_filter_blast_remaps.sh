#!/bin/bash

# Prefilter the blast alignment (basically everything except for the synteny check and needle score check) in awk. This is done to reduce the number of alignments that need to be processed by the synteny check and needle score check. The output is a file with the same format as the input file, but with only the alignments that pass the prefilter.
inputBLASTFile=$1

echo "Filtering BLAST file: $inputBLASTFile"
echo "Total number of alignments before pre-filter: $(wc -l "$inputBLASTFile" | awk '{print $1}')"

# Pre-Filter 1: Remove all entries where the first base pair of the subject sequence is not the same as the last base pair of the query sequence (post global alignment) (column 34)
# Output all rows without 0 quality score
echo "Pre-Filter 1: Remove all entries where the first base pair of the subject sequence is not the same as the last base pair of the query sequence"
awk '{if ($34 == 1) print $0}' "$inputBLASTFile" >blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_1_applied.tsv
echo "Pre-Filter 1: Done"
echo "Total number of alignments after pre-filter 1: $(wc -l blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_1_applied.tsv | awk '{print $1}')"

# Pre-Filter 2: Remove all entries where the last base pair of the subject sequence is not the same as the last base pair of the query sequence (post global alignment) (column 35)
echo "Pre-Filter 2: Remove all entries where the last base pair of the subject sequence is not the same as the last base pair of the query sequence"
awk '{if ($35 == 1) print $0}' blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_1_applied.tsv >blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_2_applied.tsv
echo "Pre-Filter 2 Done"
echo "Total number of alignments after pre-filter 2: $(wc -l blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_2_applied.tsv | awk '{print $1}')"

# Pre-Filter 3: Remove all entries where the first 3 base pairs of the subject sequence is not the same as the first base pair of the query sequence if the sequence is a CDS or exon (post global alignment) (use 3 instead of one to account for indels that could have occurred and caused frameshifts in the CDS) (column 30 for first 3 alignment) (column 42 for feature type)
echo "Pre-Filter 3: Remove all entries where the first 3 base pairs of the subject sequence is not the same as the first base pair of the query sequence if the sequence is a CDS or exon"
awk '{if (($30 == 3) && ($42 == "CDS" || $42 == "exon")) print $0; else if (! ($42 == "CDS" || $42 == "exon"))  print $0}' blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_2_applied.tsv >blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_3_applied.tsv
echo "Pre-Filter 3 Done"
echo "Total number of alignments after pre-filter 3: $(wc -l blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_3_applied.tsv | awk '{print $1}')"

# Pre-Filter 4: Remove all entries where the last 3 base pairs of the subject sequence is not the same as the last base pair of the query sequence if the sequence is a CDS or exon (post global alignment) (column 31 for last 3 alignment) (column 42 for feature type)
echo "Pre-Filter 4: Remove all entries where the last 3 base pairs of the subject sequence is not the same as the last base pair of the query sequence if the sequence is a CDS or exon"
awk '{if (($31 == 3) && ($42 == "CDS" || $42 == "exon")) print $0; else if (! ($42 == "CDS" || $42 == "exon"))  print $0}' blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_3_applied.tsv >blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_4_applied.tsv
echo "Pre-Filter 4 Done"
echo "Total number of alignments after pre-filter 4: $(wc -l blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_4_applied.tsv | awk '{print $1}')"

# Pre-Filter 5: Remove all entries where the first 8 base pairs of the subject sequence are not the same as the first 8 base pairs of the query sequence (post global alignment) (column 16)
echo "Pre-Filter 5: Remove all entries where the first 8/10 base pairs of the subject sequence are not the same as the first 8 base pairs of the query sequence"
awk '{if ($16 >= 8) print $0}' blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_4_applied.tsv >blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_5_applied.tsv
echo "Pre-Filter 5 Done"
echo "Total number of alignments after pre-filter 5: $(wc -l blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_5_applied.tsv | awk '{print $1}')"

# Pre-Filter 6: Remove all entries where the last 8/10 base pairs of the subject sequence are not the same as the last 8 base pairs of the query sequence (post global alignment) (column 17)
echo "Pre-Filter 6: Remove all entries where the last 8/10 base pairs of the subject sequence are not the same as the last 8 base pairs of the query sequence"
awk '{if ($17 >= 8) print $0}' blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_5_applied.tsv >blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_6_applied.tsv
echo "Pre-Filter 6 Done"
echo "Total number of alignments after pre-filter 6: $(wc -l blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_6_applied.tsv | awk '{print $1}')"

# Pre-Filter 7: Remove all entries where the length of the subject sequence is 50% smaller than the length of the query sequence (query length is column 36, subject length is column 50)
echo "Pre-Filter 7: Remove all entries where the length of the subject sequence is smaller than the length of the query sequence"
awk '{if ($54 >= $36 / 1.5) print $0}' blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_6_applied.tsv >blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_7_applied.tsv
echo "Pre-Filter 7 Done"
echo "Total number of alignments after pre-filter 7: $(wc -l blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_7_applied.tsv | awk '{print $1}')"

# Pre-Filter 8: Remove all entries where the length of the subject sequence is 50% larger than the length of the query sequence (query length is column 36, subject length is column 50)
echo "Pre-Filter 8: Remove all entries where the length of the subject sequence is larger than the length of the query sequence"
awk '{if ($54 <= $36 * 1.5) print $0}' blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_7_applied.tsv >blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_8_applied.tsv
echo "Pre-Filter 8 Done"
echo "Total number of alignments after pre-filter 8: $(wc -l blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_8_applied.tsv | awk '{print $1}')"

# Pre-Filter 9: Remove all entries where the length of the subject sequence is 10% smaller than the length of the query sequence if the query is a CDS or exon (query length is column 36, subject length is column 50, feature column is column 42)
echo "Pre-Filter 9: Remove all entries where the length of the subject sequence is 10% smaller than the length of the query sequence if the query is a CDS or exon"
awk '{if (($54 >= $36 / 1.1) && ($42 == "CDS" || $42 == "exon")) print $0; else if (! ($42 == "CDS" || $42 == "exon"))  print $0}' blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_8_applied.tsv >blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_9_applied.tsv
echo "Pre-Filter 9 Done"
echo "Total number of alignments after pre-filter 9: $(wc -l blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_9_applied.tsv | awk '{print $1}')"

# Pre-Filter 10: Remove all entries where the length of the subject sequence is 10% larger than the length of the query sequence if the query is a CDS or exon (query length is column 36, subject length is column 50, feature column is column 42)
echo "Pre-Filter 10: Remove all entries where the length of the subject sequence is 10% larger than the length of the query sequence if the query is a CDS or exon"
awk '{if (($54 <= $36 * 1.1) && ($42 == "CDS" || $42 == "exon")) print $0; else if (! ($42 == "CDS" || $42 == "exon"))  print $0}' blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_9_applied.tsv >blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_10_applied.tsv
echo "Pre-Filter 10 Done"
echo "Total number of alignments after pre-filter 10: $(wc -l blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_10_applied.tsv | awk '{print $1}')"

# Pre-filter 12: Filter out features that mapped in the wrong orientation
## Column 7 is the original strand orientation
## Column 60 is how the chromsome the feature derived from mapped to the new genome
## Column 22 is the new strand orientation
## If the original strand is +

cp blast_outputs/intermediate_filters/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_filter_10_applied.tsv ./blast_outputs/all_blasted_to_new_assembly_WITH_ALL_INFO_FINAL_and_PRE_FILTERED.tsv

# rm blast_outputs/intermediate_filters/*

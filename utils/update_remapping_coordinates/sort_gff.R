#!/usr/bin/env Rscript

# Load Required Libraries
suppressWarnings(library(tidyverse, quietly = TRUE))

# Import GFF file
commandLineArgsInputFilePath <- commandArgs(trailingOnly = TRUE)[1]
print(paste("Importing GFF...", commandLineArgsInputFilePath))
gff <- read_tsv(commandLineArgsInputFilePath, col_names = FALSE)
colnames(gff) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "q_metadata")


# Sort Dataframe
# Sort the dataframe by the first column (chromosome) (numerically after the "_") and then by the fourth column (start position)
# 1. First rename the seqid column to everything after the "_" and convert to numeric
# 2. Then sort the dataframe by the new seqid column and the start column
# 3. Then rename the seqid column back to the original name (convert to character and add the "chromosome_" back on)

print("Sorting GFF")
gff <- gff %>% mutate(seqid = as.numeric(str_sub(seqid, str_locate(seqid, "_")[2,2] + 1, -1)))
gff <- gff %>% arrange(seqid, start)
gff <- gff %>% mutate(seqid = as.character(paste0("chromosome_", seqid)))
print("Sorted GFF")

# Write to File
# Write the sorted dataframe to a new file

outFile <- "genome_data/new/final_remapped_data_updated_coordinates_sorted_minus_cp_and_mt_minus_original_metadata.gff"
print(paste("Writing to file...", outFile))
write_tsv(gff, outFile, col_names = FALSE)


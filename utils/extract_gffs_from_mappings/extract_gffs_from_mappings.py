# Import necessary modules
import sys
import argparse
import numpy as np
import pandas as pd
import re

# Just takes in all the arguments and sets up the help text
parser = argparse.ArgumentParser(description='Filters a PAF file')
parser.add_argument("-p", "--pafFile", help = "Input PAF file")
parser.add_argument("-g", "--gffFile", help = "Input GFF file of old features")
parser.add_argument("-o", "--outputFile", help = "Output file")


args = parser.parse_args()

try:
    paf = args.pafFile
    gff = args.gffFile
    outFile = args.outputFile
except:
    print("Please provide all arguments")
    sys.exit(1)
    
separator = "\t"

# import paf file as a pandas dataframe
paf_df = pd.read_csv(paf, sep = separator, header = None)

# Rename the column names of the paf file 
paf_df.columns = ["qname", "qlen", "qstart", "qend", "strand", "tname", "tlen", "tstart", "tend", "nmatch", "alen", "mapq"]

# add id column to paf file
paf_df['id'] = paf_df.index

# import gff file as a pandas dataframe
gff_df = pd.read_csv(gff, sep = separator, header = None)

# Rename the columns of the GFF file
gff_df.columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes", "featureLength"]


# loop through the gff file 
for indexi, rowi in gff_df.iterrows():
    # extract the chromosome name
    chrom = str(rowi[0])
    # extract the strand
    strand = rowi[6]
    # If the strand is negative, the start and end positions need to be swapped
    if strand == "-":
        # extract the start position and end position for the negative strand
        start = rowi[4]
        end = rowi[3]
    elif strand == "+":
        # extract the start position and end position for the positive strand
        start = rowi[3]
        end = rowi[4]
    else: 
        gff_df.loc[indexi, 'id'] = -1
        gff_df.loc[indexi, 'target_chrom'] = -1
        gff_df.loc[indexi, 'mapping_start_on_old'] = -1
        gff_df.loc[indexi, 'mapping_end_on_old'] = -1
        gff_df.loc[indexi, 'if_query_and_target_on_same_strand'] = -1 # ‘+’ if query/target on the same strand; ‘-’ if opposite
        gff_df.loc[indexi, 'mapping_start_on_new'] = -1
        gff_df.loc[indexi, 'mapping_end_on_new'] = -1
        gff_df.loc[indexi, 'mapping_q_score'] = -1
        continue

    # Loop through the paf file
    for indexj, rowj in paf_df.iterrows():
        # check if the chromosome name matches
        if rowj[0] == chrom:
            # check if the start and end positions from the gff file are within the start and end positions from the paf file
            # the paf file is 0 based, so the start position needs to be subtracted by 1
            if start - 1 >= rowj[2] and end <= rowj[3]:
                # add the id column to the gff file
                gff_df.loc[indexi, 'id'] = rowj['id']
                # add the target chromosome name to the gff file
                gff_df.loc[indexi, 'target_chrom'] = rowj[5]
                # add the regions from the original chromosome that the feature belongs to 
                gff_df.loc[indexi, 'mapping_start_on_old'] = rowj[2]
                gff_df.loc[indexi, 'mapping_end_on_old'] = rowj[3]
                # add the strand information
                gff_df.loc[indexi, 'if_query_and_target_on_same_strand'] = rowj[4] # ‘+’ if query/target on the same strand; ‘-’ if opposite
                # add the regions from the new chromosome that the feature should belong to
                gff_df.loc[indexi, 'mapping_start_on_new'] = rowj[7]
                gff_df.loc[indexi, 'mapping_end_on_new'] = rowj[8]
                # add the mapping quality score of the old assembly to new assembly alignment
                gff_df.loc[indexi, 'mapping_q_score'] = rowj[11]
                break
    # If no match is found, set the id column to -1
    if not gff_df.loc[indexi, 'id'] >= 0:
        gff_df.loc[indexi, 'id'] = -1
        gff_df.loc[indexi, 'target_chrom'] = -1
        gff_df.loc[indexi, 'mapping_start_on_old'] = -1
        gff_df.loc[indexi, 'mapping_end_on_old'] = -1
        gff_df.loc[indexi, 'if_query_and_target_on_same_strand'] = -1 
        gff_df.loc[indexi, 'mapping_start_on_new'] = -1
        gff_df.loc[indexi, 'mapping_end_on_new'] = -1
        gff_df.loc[indexi, 'mapping_q_score'] = -1

    if indexi % 100 == 0 or indexi == len(gff_df):
        print(f"Row {indexi} of {len(gff_df)} complete")

# Remove decimals from added columns
gff_df['id'] = gff_df['id'].astype(int)
gff_df['target_chrom'] = gff_df['target_chrom'].astype(str).apply(lambda x: re.sub( r'\.0$', '', x) )
gff_df['mapping_start_on_old'] = gff_df['mapping_start_on_old'].astype(int)
gff_df['mapping_end_on_old'] = gff_df['mapping_end_on_old'].astype(int)
gff_df['mapping_start_on_new'] = gff_df['mapping_start_on_new'].astype(int)
gff_df['mapping_end_on_new'] = gff_df['mapping_end_on_new'].astype(int)
gff_df['mapping_q_score'] = gff_df['mapping_q_score'].astype(int)

# Output the gff file
print(f"Writing to file...{outFile}")
pd.DataFrame.to_csv(gff_df, outFile, sep = separator, header = False, index = False)
            
        
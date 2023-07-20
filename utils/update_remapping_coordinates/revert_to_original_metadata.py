# Import required modules
import sys
import argparse
import numpy as np
import pandas as pd
import pickle as pkl

# Import data
parser = argparse.ArgumentParser(description='Replaces the metadata in a GFF file with the orignial metadata')
parser.add_argument("-g", "--gffFile", help = "GFF file")
parser.add_argument("-m", "--metaTable", help = "Hash table of the new metadata to the old metadata")
parser.add_argument("-o", "--outputGffFile", help = "Output file with updated metadata")


args = parser.parse_args()

separator = "\t"

try:
    gffFile = args.gffFile
    metaTable = args.metaTable
    outputGffFile = args.outputGffFile

except:
    print("Please provide all arguments")
    sys.exit(1)

gff = pd.read_csv(gffFile, sep = '\t', header = None)
metadata_table = pd.read_csv(metaTable, sep = '\t', header = None)

# Create a dictionary (hash table) from the metadata table
print("Creating hash table from metadata")
metadata_dict = {}
counter = 0
for row in metadata_table.iterrows():
  metadata_dict[row[1][0]] = row[1][1]
  
  # Print progress
  counter += 1
  if counter % 1000 == 0:
    print(f'Processed {counter} of {len(metadata_table)}')
  
  
# Replace the metadata in the GFF file
print("Replacing cleaned up metadata with original")
counter = 0
for row in gff.iterrows():
  if row[1][8] in metadata_dict:
    # replace the metadata in the original table
    gff.iloc[row[0], 8] = metadata_dict[row[1][8]]
    
  # Print progress
  counter += 1
  if counter % 1000 == 0:
    print(f'Processed {counter} of {len(gff)}')
    
# Write the new GFF file
print(f"Writing new GFF file to {outputGffFile}")
pd.DataFrame(gff).to_csv(outputGffFile, sep = separator, header = False, index = False)

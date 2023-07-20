# Import required modules
import sys
import argparse
import numpy as np
import pandas as pd
import pickle as pkl

# Import data
parser = argparse.ArgumentParser(description='Filters a BLAST file post global alignment for proper remaps')
# parser.add_argument("-g", "--gffFile", help = "Input Original GFF file")
parser.add_argument("-b", "--blastFile", help = "Input BLAST file post global alignment and all other checks")
parser.add_argument("-p", "--pafFile", help = "Input PAF of old to new assembly alignment with no gaps")
parser.add_argument("-o", "--outFile", help = "Output file for blast table filtered post synteny check")
# parser.add_argument("-or", "--outputRemappedFile", help = "Output file for remapped features")
# parser.add_argument("-on", "--outputNonRemappedFile", help = "Output file for non remapped features")
parser.add_argument("-ns", "--noSave", action = "store_true", help = "Don't save pkl files (default is to save)")
parser.add_argument("-nl", "--noLoad", action = "store_true", help = "Don't load pkl files (default is to load)")

args = parser.parse_args()

separator = "\t"

try:
    # oldGffFile = args.gffFile
    blastFile = args.blastFile
    oldToNewAlignmentPafFile = args.pafFile
    outFile = args.outFile
    # outputRemappedFile = args.outputRemappedFile
    # outputNonRemappedFile = args.outputNonRemappedFile
    
except:
    print("Please provide all arguments")
    sys.exit(1)

# oldGff = pd.read_csv(oldGffFile, sep = '\t', header = None)
postGlobalAlignmentBlast = pd.read_csv(blastFile, sep = '\t', header = None)
oldToNewAlignmentPaf = pd.read_csv(oldToNewAlignmentPafFile, sep = '\t', header = None)

# Strip last empty column from oldGff
# oldGff = oldGff.iloc[:, :-1]
# Rename the columns of the GFF file
# oldGff.columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

# Rename the column names of the BLAST file
postGlobalAlignmentBlast.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sstrand", "order", "needle_score", "first_10", "last_10", "first_9", "last_9", "first_8", "last_8", "first_7", "last_7", "first_6", "last_6", "first_5", "last_5", "first_4", "last_4", "first_3", "last_3", "first_2", "last_2", "first_1", "last_1", "q_length", "q_chr", "q_old_start", "q_old_end", "q_old_strand", "q_old_source", "q_old_type", "q_old_score", "q_old_phase", "q_metadata", "mapping_row_in_paf_0i", "expected_chromosome", "mapping_start_on_old_chromosome", "mapping_end_on_old_chromosome", "if_query_and_target_on_same_strand", "mapping_start_on_new_chromosome", "mapping_end_on_new_chromosome", "mapping_q_score", "s_length"]


# Rename the column names of the oldToNewAlignmentPaf 
oldToNewAlignmentPaf.columns = ["qname", "qlen", "qstart", "qend", "strand", "tname", "tlen", "tstart", "tend", "nmatch", "alen", "mapq"]

# Grab the unique query sequence IDs from the BLAST file
uniqueQuerySequenceIDs = postGlobalAlignmentBlast["qseqid"].unique()

# # function to make acceptable region list
# def make_acceptable_region_dict():
    
#     # Make a vector of character names
#     new_chromosomes = np.char.add("chromosome_", np.arange(1, 26).astype(str))
#     new_chromosomes = np.append(new_chromosomes, ["mitochondrion","chloroplast"])
    
#     # Make an empty vector of NA's
#     tmp = np.repeat(np.nan, 27)
    
#     # Make a dictionary from the 2 vectors (make 2 lists, one adhering to heuristic and one adhering to judgment calls)
#     acceptable_regions = dict(zip(new_chromosomes, tmp))
    
#     # add chromosome, start, and end to each entry of the list
#     labels = ["chromosome", "start", "end", "old_start", "old_end","strand"]
#     for chrom in new_chromosomes:
#         acceptable_regions[chrom] = dict(zip(labels, np.repeat(np.nan, len(labels))))
    
#     # return the final dictionary
#     return acceptable_regions

# # function to replace NA's with acceptable regions
# def add_acceptable_region(acc_reg_dict, new, old, start, end, old_start, old_end, strand):
    
#     # # if this is the first acceptable region added
#     try:
#       if np.isnan(acc_reg_dict[new]["chromosome"]):
#           acc_reg_dict[new]["chromosome"] = str(old)
#           acc_reg_dict[new]["start"] = start
#           acc_reg_dict[new]["end"] = end
#           acc_reg_dict[new]["old_start"] = old_start
#           acc_reg_dict[new]["old_end"] = old_end
#           acc_reg_dict[new]["strand"] = strand
#           return acc_reg_dict
#     except TypeError:
#         # if this is not the first acceptable region added
#         acc_reg_dict[new]["chromosome"] = np.append(acc_reg_dict[new]["chromosome"], str(old))
#         acc_reg_dict[new]["start"] = np.append(acc_reg_dict[new]["start"], start)
#         acc_reg_dict[new]["end"] = np.append(acc_reg_dict[new]["end"], end)
#         acc_reg_dict[new]["old_start"] = np.append(acc_reg_dict[new]["old_start"], old_start)
#         acc_reg_dict[new]["old_end"] = np.append(acc_reg_dict[new]["old_end"], old_end)
#         acc_reg_dict[new]["strand"] = np.append(acc_reg_dict[new]["strand"], strand)
    
#     return acc_reg_dict

# # initialize lists
# acceptable_regions = make_acceptable_region_dict()

# # Run the add_acceptable_region function for every row of the oldToNewAlignmentPaf file
# for index, row in oldToNewAlignmentPaf.iterrows():
#     acceptable_regions = add_acceptable_region(acc_reg_dict=acceptable_regions, new=row["tname"], old=row["qname"], start=int(row["tstart"]), end=int(row["tend"]), old_start=row["qstart"], old_end=row["qend"], strand=row["strand"])

# manually add acceptable regions for mitochondrion and chloroplast

## mitochondrion
# # Pt mitochondria is 77356 bp long
# #TODO check if "mitochondria" or "mitochondrion" is correct
# acceptable_regions = add_acceptable_region(acc_reg_dict=acceptable_regions, new="mitochondrion", old="mitochondria", start=1, end=77356, old_start=1, old_end=77356, strand = "+")

# #chloroplast
# # Pt chloroplast is 117369 bp long
# acceptable_regions = add_acceptable_region(acc_reg_dict=acceptable_regions, new="chloroplast", old="chloroplast", start=1, end=117369, old_start=1, old_end=117369, strand = "+")

# Apply filters to the BLAST file

# Filter 1: Remove all entries where the first base pair of the subject sequence is not the same as the last base pair of the query sequence (post global alignment)
# postGlobalAlignmentBlast = postGlobalAlignmentBlast[postGlobalAlignmentBlast["first_3"] == 3]

# Filter 2: Remove all entries where the last base pair of the subject sequence is not the same as the last base pair of the query sequence (post global alignment)


# Filter 5: Remove all entries where the first 8 base pairs of the subject sequence are not the same as the first 8 base pairs of the query sequence (post global alignment)
# bp_cutoff = 8
# postGlobalAlignmentBlast = postGlobalAlignmentBlast[postGlobalAlignmentBlast["first_10"] >= bp_cutoff]

# Filter 6: Remove all entries where the last 8/10 base pairs of the subject sequence are not the same as the last 8 base pairs of the query sequence (post global alignment)
# postGlobalAlignmentBlast = postGlobalAlignmentBlast[postGlobalAlignmentBlast["last_10"] >= bp_cutoff]

# Filter 7: Remove all entries where the length of the subject sequence is smaller than the length of the query sequence 
# length_cutoff = 1.5
# postGlobalAlignmentBlast = postGlobalAlignmentBlast[postGlobalAlignmentBlast["s_length"] >= postGlobalAlignmentBlast["q_length"] / length_cutoff]

# Filter 8: Remove all entries where the length of the subject sequence is larger than the length of the query sequence
# postGlobalAlignmentBlast = postGlobalAlignmentBlast[postGlobalAlignmentBlast["s_length"] <= postGlobalAlignmentBlast["q_length"] * length_cutoff]

# Filter 9: Remove all entries where the subject sequence is not in the acceptable region list
# Create a function to loop through the BLAST file and add a column to the BLAST file that indicates whether the subject sequence is in the acceptable region list
# def syneny_check(BLASTtable, acceptable_regions):
#   # Add column to end of BLAST table with empty values
#     BLASTtable["in_acceptable_region"] = np.repeat(np.nan, len(BLASTtable))
#     BLASTtable["how_features_old_chr_maps_to_new"] = np.repeat(np.nan, len(BLASTtable))
#     BLASTtable["orientation_check"] = np.repeat(np.nan, len(BLASTtable))
    
    
#     loop_count = 0
#     # loop through each row of the blast table
#     for indexi, row in BLASTtable.iterrows():
#       # extract the chromosome that the BLAST entry mapped to
#       s_chr = row["sseqid"]
#       # Extract the strand of the BLAST entry
#       s_strand = row["sstrand"]
#       # Extract the start position of the BLAST entry
#       s_start = int(row["sstart"])
#       # Extract the end position of the BLAST entry
#       s_end = int(row["send"])
#       # Extract the original chromosome that the BLAST entry mapped to
#       q_chr = str(row["q_chr"])
#       # Extract the old start position of the BLAST entry
#       q_start = int(row["q_old_start"])
#       # Extract the old end position of the BLAST entry
#       q_end = int(row["q_old_end"])
#       # Extract the old strand of the BLAST entry
#       q_strand = row["q_old_strand"]
      
      
#       # Loop through the list of acceptable regions for the subject chromosome
#       for indexj, old_chr in enumerate(acceptable_regions[s_chr]["chromosome"]):
        
#         # If the subject chromosome has more than one mapping, then extract the start and end positions of the acceptable region from an array. Otherwise extract the start and end positions from a single value.
#         if isinstance(acceptable_regions[s_chr]["chromosome"], np.ndarray):
#           acceptable_start = acceptable_regions[s_chr]["start"][indexj]
#           acceptable_end = acceptable_regions[s_chr]["end"][indexj]
#           acceptable_old_start = acceptable_regions[s_chr]["old_start"][indexj]
#           acceptable_old_end = acceptable_regions[s_chr]["old_end"][indexj]
#           old_to_new_mapping = acceptable_regions[s_chr]["strand"][indexj]
#         else:
#           acceptable_start = acceptable_regions[s_chr]["start"]
#           acceptable_end = acceptable_regions[s_chr]["end"]
#           acceptable_old_start = acceptable_regions[s_chr]["old_start"]
#           acceptable_old_end = acceptable_regions[s_chr]["old_end"]
#           old_to_new_mapping = acceptable_regions[s_chr]["strand"]
          
#         # Check if the subject chromosome is in the acceptable region list
#         if old_chr == q_chr:
#           # check if the start position of the BLAST entry is greater than the start position of the acceptable region
#           if s_start >= acceptable_start:
#             # check if the end position of the BLAST entry is less than the end position of the acceptable region
#             if s_end <= acceptable_end:
#               #  check if the start position of the BLAST entry query is greater than the start position of the region of the old chromosome mapping that the query sequence should fall within
#               if q_start >= acceptable_old_start:
#                 # check if the end position of the BLAST entry query is less than the end position of the region of the old chromosome mapping that the query sequence should fall within
#                 if q_end <= acceptable_old_end:
#                   # if all of the above conditions are met, then the BLAST entry is in the acceptable region list
#                   BLASTtable.loc[indexi, "in_acceptable_region"] = True
#                   BLASTtable.loc[indexi, "how_features_old_chr_maps_to_new"] = old_to_new_mapping
                  
#                   # Now check if mapped in the correct orientation
#                   ## Column 7 is the original strand orientation
#                   ## Column 60 is how the chromsome the feature derived from mapped to the new genome
#                   ## Column 22 is the new strand orientation
#                   ## If the original strand is + and the original chromosome mapped to the new genome is +, then the new strand should be +
#                   if q_strand == "+" and old_to_new_mapping == "+":
#                     if s_strand == "plus":
#                       BLASTtable.loc[indexi, "orientation_check"] = True
#                     elif s_strand == "minus":
#                       BLASTtable.loc[indexi, "orientation_check"] = False
#                   ## If the original strand is - and the original chromosome mapped to the new genome is +, then the new strand should be -
#                   if q_strand == "-" and old_to_new_mapping == "+":
#                     if s_strand == "plus":
#                       BLASTtable.loc[indexi, "orientation_check"] = False
#                     elif s_strand == "minus":
#                       BLASTtable.loc[indexi, "orientation_check"] = True
#                   ## If the original strand is + and the original chromosome mapped to the new genome is -, then the new strand should be -
#                   if q_strand == "+" and old_to_new_mapping == "-":
#                     if s_strand == "plus":
#                       BLASTtable.loc[indexi, "orientation_check"] = False
#                     elif s_strand == "minus":
#                       BLASTtable.loc[indexi, "orientation_check"] = True
#                   ## If the original strand is - and the original chromosome mapped to the new genome is -, then the new strand should be +
#                   if q_strand == "-" and old_to_new_mapping == "-":
#                     if s_strand == "plus":
#                       BLASTtable.loc[indexi, "orientation_check"] = True
#                     elif s_strand == "minus":
#                       BLASTtable.loc[indexi, "orientation_check"] = False
                      
#                   # Break the loop
#                   break
            
#       # Otherwise, the BLAST entry is not in the acceptable region list. If the BLAST entry has not been marked true by this point, set it to false
#       if pd.isnull(BLASTtable.loc[indexi, "in_acceptable_region"]):
#         BLASTtable.loc[indexi, "in_acceptable_region"] = False
        
#       # Print the progress of the loop
#       loop_count += 1
#       if loop_count % 1000 == 0:
#         print(f"Finished {loop_count} of {len(BLASTtable.index)}")
      
#     # Return the BLAST table with the new columns after the loop has finished
#     return BLASTtable


# Modified function based on new blast table that has the assembly mapping information
def syneny_check(BLASTtable):
  # Add column to end of BLAST table with empty values
    BLASTtable["in_acceptable_region"] = np.repeat(np.nan, len(BLASTtable))
    BLASTtable["orientation_check"] = np.repeat(np.nan, len(BLASTtable))
    
    
    loop_count = 0
    # loop through each row of the blast table
    for indexi, row in BLASTtable.iterrows():
      # extract the chromosome that the BLAST entry mapped to
      s_chr = row["sseqid"]
      
      # Extract the strand of the BLAST entry
      ## Either "plus" or "minus"
      s_strand = row["sstrand"]
      ## Convert "plus" to "+" and "minus" to "-"
      if s_strand == "plus":
        s_strand = "+"
      elif s_strand == "minus":
        s_strand = "-"
        
      # Extract the start position of the BLAST entry
      s_start = int(row["sstart"])
      # Extract the end position of the BLAST entry
      s_end = int(row["send"])
      # Extract the original chromosome that the BLAST entry mapped to
      q_chr = str(row["q_chr"])
      # Extract the old start position of the BLAST entry
      q_start = int(row["q_old_start"])
      # Extract the old end position of the BLAST entry
      q_end = int(row["q_old_end"])
      # Extract the old strand of the BLAST entry
      ## Either "+" or "-"
      q_strand = row["q_old_strand"]
      # Extract the expected mapping of the old chromosome to the new chromosome
      expected_chromosome = row["expected_chromosome"]
      # Extract the start position of the old chromosome mapping to the old assembly chromosome or scaffold
      mapping_start_on_old_chromosome = int(row["mapping_start_on_old_chromosome"])
      # Extract the end position of the old chromosome mapping to the old assembly chromosome or scaffold
      mapping_end_on_old_chromosome = int(row["mapping_end_on_old_chromosome"])
      # Extract whether the old chromosome mapping and new chromosome are on the same strand
      ## "+" if they are on the same strand, "-" if they are on the opposite strand
      if_query_and_target_on_same_strand = row["if_query_and_target_on_same_strand"]
      # Extract the start position of the old chromosome mapping to the new assembly chromosome
      mapping_start_on_new_chromosome = int(row["mapping_start_on_new_chromosome"])
      # Extract the end position of the old chromosome mapping to the new assembly chromosome
      mapping_end_on_new_chromosome = int(row["mapping_end_on_new_chromosome"])
          
        # Check if the subject chromosome is where the BLAST entry should map to
      if s_chr == expected_chromosome:
        # check if the start position of the BLAST entry is greater than the start position of the acceptable region
        if s_start >= mapping_start_on_new_chromosome:
          # check if the end position of the BLAST entry is less than the end position of the acceptable region
          if s_end <= mapping_end_on_new_chromosome:
            #  check if the start position of the BLAST entry query is greater than the start position of the region of the old chromosome mapping that the query sequence should fall within
            if q_start >= mapping_start_on_old_chromosome:
              # check if the end position of the BLAST entry query is less than the end position of the region of the old chromosome mapping that the query sequence should fall within
              if q_end <= mapping_end_on_old_chromosome:
                # if all of the above conditions are met, then the BLAST entry is in the acceptable region list
                BLASTtable.loc[indexi, "in_acceptable_region"] = True
                
                # Now check if the BLAST entry is mapped to the correct orientation
                # if_query_and_target_on_same_strand was extracted from the old assembly to new assembly mapping .paf file
                ## in .paf files, "+" if they are on the same strand, "-" if they are on the opposite strand
                
                # If the old chromosome mapping and new chromosome are on the same strand, then the BLAST entry should be on the same as the query sequence
                if if_query_and_target_on_same_strand == "+":
                  if q_strand == s_strand:
                    BLASTtable.loc[indexi, "orientation_check"] = True
                  elif q_strand != s_strand:
                    BLASTtable.loc[indexi, "orientation_check"] = False
                
                # If the old chromosome mapping and new chromosome are on the opposite strand, then the BLAST entry should be on the opposite strand as the query sequence
                elif if_query_and_target_on_same_strand == "-":
                  if q_strand == s_strand:
                    BLASTtable.loc[indexi, "orientation_check"] = False
                  elif q_strand != s_strand:
                    BLASTtable.loc[indexi, "orientation_check"] = True


      # Otherwise, the BLAST entry is not in the acceptable region list. If the BLAST entry has not been marked true by this point, set it to false
      if pd.isnull(BLASTtable.loc[indexi, "in_acceptable_region"]):
        BLASTtable.loc[indexi, "in_acceptable_region"] = False
        BLASTtable.loc[indexi, "orientation_check"] = False
        
      # Print the progress of the loop
      loop_count += 1
      if loop_count % 1000 == 0:
        print(f"Finished {loop_count} of {len(BLASTtable.index)}")
      
    # Return the BLAST table with the new columns after the loop has finished
    return BLASTtable


if args.noLoad:
  postGlobalAlignmentBlastPostSyntenyCheck = syneny_check(postGlobalAlignmentBlast)
else:
    try:
      filePath = "./blast_outputs/BLAST_post_synteny_check.pkl"
      with open(filePath, 'rb') as f:
        print(f"Loading postGlobalAlignmentBlastPostSyntenyCheck from {filePath}...")
        postGlobalAlignmentBlastPostSyntenyCheck = pkl.load(f)
        print(f"Loaded")
    except:
        print(f"Loading {filePath} failed. Performing synteny check on postGlobalAlignmentBlast from scratch...")
        postGlobalAlignmentBlastPostSyntenyCheck = syneny_check(postGlobalAlignmentBlast)
        if not args.noSave:
          filePath = "./blast_outputs/BLAST_post_synteny_check.pkl"
          with open(filePath, 'wb') as f:
            print(f"Saving postGlobalAlignmentBlastPostSyntenyCheck to {filePath}...")
            pkl.dump(postGlobalAlignmentBlastPostSyntenyCheck, f)
            print(f"Saved")


# Perform filtering
# Filter 1: Filter out features that mapped to the wrong chromosome
print("Filtering out features that mapped to the wrong chromosome...")
print(f"The number of features before filtering is {len(postGlobalAlignmentBlastPostSyntenyCheck.index)}")
postGlobalAlignmentBlastPostSyntenyCheck = postGlobalAlignmentBlastPostSyntenyCheck[postGlobalAlignmentBlastPostSyntenyCheck["in_acceptable_region"] == True]
print("Filtering 1 complete")
print(f"The number of features after filtering is {len(postGlobalAlignmentBlastPostSyntenyCheck.index)}")

# Filter 2: Filter out features that mapped in the wrong orientation
print("Filtering out features that mapped in the wrong orientation...")
print(f"The number of features before filtering is {len(postGlobalAlignmentBlastPostSyntenyCheck.index)}")
postGlobalAlignmentBlastPostSyntenyCheck = postGlobalAlignmentBlastPostSyntenyCheck[postGlobalAlignmentBlastPostSyntenyCheck["orientation_check"] == True]
print("Filtering 2 complete")
print(f"The number of features after filtering is {len(postGlobalAlignmentBlastPostSyntenyCheck.index)}")


# Save the BLAST table after filtering
if not args.noSave:
  print(f"Saving BLAST table after filtering to {outFile}...")
  postGlobalAlignmentBlastPostSyntenyCheck.to_csv(outFile, sep="\t", index=False)

# Filter 8: Filter by highest needle score
# Perform this in R


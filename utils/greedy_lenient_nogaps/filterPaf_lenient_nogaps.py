import argparse
import sys
import pickle as pkl

# Just takes in all the arguments and sets up the help text
parser = argparse.ArgumentParser(description='Filters a PAF file')
parser.add_argument("-f1", "--file1", help = "Input paf file")
parser.add_argument("-gs", "--genomeSizes", help = "Input chromosome size file")
parser.add_argument("-of", "--outputFile", help = "Output file")
parser.add_argument("-s", "--sensitivity", help = "Input a sensitivity for overlaps with other fragements (ranging from 0-1)")
parser.add_argument("-ns", "--noSave", action = "store_true", help = "Don't save pkl files (default is to save)")
parser.add_argument("-nl", "--noLoad", action = "store_true", help = "Don't load pkl files (default is to load)")

args = parser.parse_args()

separator = "\t"

try:
    file1 = args.file1
    genSizes = args.genomeSizes
    outFile = args.outputFile
    acceptableOverlap = float(args.sensitivity)
except:
    print("Please provide all arguments")
    sys.exit(1)
    
# Read in genome size + put in dict
names = []
size = []
with open(genSizes,'r') as f:
    lines = f.readlines();
    for line in lines:
        lineVal = line.rstrip()
        parts = lineVal.split(separator)
        names.append(parts[0].rstrip().split("_")[-1])
        size.append(parts[1].rstrip())

mappingInfoOnChromosomes = {k: [] for k in names}

def buildChromosomes(mappingInfoOnChromosomes, size):
    # Put the ranges inside the dict with the keys being the names (1,2,3....chloroplast,mitochondrion)
    # 1's indicate this position hasnt been "used up" yet
    # 0's indicate that it has already been matched to
    for index, key in enumerate(mappingInfoOnChromosomes.keys()):
        for x in range(0,int(size[index])):
            mappingInfoOnChromosomes[key].append({"available": True, "overlap": False, "target_position": x, "chromosomes_matched": []})
            if x % 100000 == 0:
                print(f"chromosome {key}: {x} of {size[index]} bases added to mappingInfoOnChromosomes dictionary ")
    return mappingInfoOnChromosomes

# Try to load the mappingInfoOnChromosomes dict from a pickle file if it exists to save time
if args.noLoad:
    mappingInfoOnChromosomes = buildChromosomes(mappingInfoOnChromosomes, size)
else:
    try:
        filePath = "./chr_sizes/mappingInfoOnChromosomes.pkl"
        with open(filePath, 'rb') as f:
            print(f"loading mappingInfoOnChromosomes from {filePath}...")
            mappingInfoOnChromosomes = pkl.load(f)
            print(f"Loaded")
    except:
        print(f"Loading {filePath} failed. Building mappingInfoOnChromosomes from scratch...")
        mappingInfoOnChromosomes = buildChromosomes(mappingInfoOnChromosomes, size)

if not args.noSave:
    filePath = "./chr_sizes/mappingInfoOnChromosomes.pkl"
    with open(filePath, 'wb') as f:
        print(f"Saving mappingInfoOnChromosomes to {filePath}...")
        pkl.dump(mappingInfoOnChromosomes, f)
        print(f"Saved")

def buildPafTable(file1, separator):
    pafTable = []
    with open(file1,'r') as f:
        lines = f.readlines();
        for i, line in enumerate(lines):
            lineVal = line.rstrip()
            parts = lineVal.split(separator)
            parts[10] = int(parts[10])
            pafTable.append(parts)
    return pafTable

# Try to load the pafTable list from a pickle file if it exists to save time
if args.noLoad:
    pafTable = buildPafTable(file1, separator)
else:
    try:
        filePath = "./minimap_outputs/pafTable.pkl"
        with open(filePath, 'rb') as f:
            print(f"Loading pafTable from {filePath}...")
            pafTable = pkl.load(f)
            print(f"Loaded")
    except:  
        print(f"Loading {filePath} failed. Building pafTable from scratch...")
        pafTable = buildPafTable(file1, separator)
        
if not args.noSave:
    filePath = "./minimap_outputs/pafTable.pkl"
    with open(filePath, 'wb') as f:
        print(f"Saving pafTable to {filePath}...")
        pkl.dump(pafTable, f)
        print(f"Saved")

# largest to smallest sort
print("Sorting pafTable by mapping length...")
pafTable.sort(key=lambda x:x[10], reverse=True)
print("Sorted")


def insideOfValidAreaOverlap(targetStart, targetEnd, chrPotMatchArea):
    length = abs(targetEnd - targetStart)
    overlappingBases = 0

    # don't need to +1 targetStart and targetEnd because .paf files are 0-based
    for i in range(targetStart, targetEnd):
        if chrPotMatchArea[i]["available"] == False:
            overlappingBases+=1
    
    percentOverlap = overlappingBases/length
    return percentOverlap 

def addMappingInfo(targetStart, targetEnd, queryChr, queryStart, queryEnd, queryStrand, chrPotMatchArea):
    if queryStrand == "-":
        locationOnQuery = queryEnd - 1 # -1 because .paf files are 0-based and end coords are non-inclusive
    else:
        locationOnQuery = queryStart
    # don't need to +1 targetStart and targetEnd because .paf files are 0-based
    for i in range(targetStart, targetEnd):
        chrPotMatchArea[i]["available"] = False
        # add chromosomes that match to this spot
        chrPotMatchArea[i]["chromosomes_matched"].append({"queryChr": queryChr, "queryMatchedToStrand": queryStrand, "locationOnQuery": locationOnQuery})
        if queryStrand == "-":
            locationOnQuery -= 1
        else:
            locationOnQuery += 1
        # check if there is an overlap
        if len(chrPotMatchArea[i]["chromosomes_matched"]) > 1:
            chrPotMatchArea[i]["overlap"] = True

# First loop: Go through the maps from largest to smallest and add them to the mappingInfoOnChromosomes dictionary if they don't overlap with anything more than the acceptableOverlap
print("Starting the first loop to add all the mappings with acceptable overlaps")

first_pass_table = []

for rowi, row in enumerate(pafTable):
    queryChr = row[0]
    queryStart = int(row[2])
    queryEnd = int(row[3])
    queryStrand = row[4]
    targetStart = int(row[7])
    targetEnd = int(row[8])
    matchingChromosome = row[5].split("_")[-1]
    
    chrPotMatchArea = mappingInfoOnChromosomes[matchingChromosome]
    
    percentOverlap = insideOfValidAreaOverlap(targetStart, targetEnd, chrPotMatchArea)

    if percentOverlap <= acceptableOverlap:
        addMappingInfo(targetStart, targetEnd, queryChr, queryStart, queryEnd, queryStrand, chrPotMatchArea)
        first_pass_table.append(row)
    
    # Print progress
    if rowi % 10 == 0:
        print(f"Progress: {rowi}/{len(pafTable)}", end="\r")
    
print("Done with the first loop")

# Save the table after the initial mapping to an intermediate file 
intermediateFilePath = "./filtered_alignments/greedy-lenient-nogaps/oldToNewMap-450size-20score_20pct-overlap_intermediate.paf"

# Build the intermediate table string
print("Formatting intermediate table for output...")
outputStringIntermediate = ""
for row in first_pass_table:
    for i, value in enumerate(row):
        if i!= len(row)-1:
            outputStringIntermediate += str(value) + "\t"
        else:
            outputStringIntermediate += str(value)
    outputStringIntermediate += "\n"
print("Done formatting intermediate table for output")

# Output the intermediate table
with open(intermediateFilePath, 'w') as f:
    print(f"Writing intermediate table to {intermediateFilePath}...")
    f.write(outputStringIntermediate)
    print("Done")

# Now that the initial mapping is done, we need to go through and remove any fragments that overlap with other fragments

# Each mapping can have to places where it overlaps with another set of mappings (one on the left and one on the right)
# Loop through each mapping and see what the mapping situation is
# There are all the following outcomes and what should happen:

# 1. The mapping has 2 regions of overlaps, one on the left and one on the right

# 2. The mapping has 1 overlap on the left

# 3. The mapping has 1 overlap on the right

# 4. The mapping has no overlaps
    # Leave the mapping as is

# 5. The entire mapping is overlapped
    # If for some reason a mapping was not properly filtered out and remains in the data, Remove the mapping 
    # Make sure it doesn't enter the final table
    
def removeOverlap(targetStart, targetEnd, queryChr, queryStart, queryEnd, chrPotMatchArea):
    newTargetStart = targetStart
    newTargetEnd = targetEnd
    newQueryStart = queryStart
    newQueryEnd = queryEnd


    # Run this loop to determine what kind of overlap the mapping has
    twoOverlapPattern = [True, False, True]
    leftOverlapPattern = [True, False]
    rightOverlapPattern = [False, True]
    noOverlapPattern = [False]
    completeOverlapPattern = [True]

    mappingOverlapPattern = {"overlapPattern": [], "overlap_region_1": {"start": 0, "end": 0}, "overlap_region_2": {"start": 0, "end": 0}}
    for i in range(targetStart, targetEnd):
        # add the information of the first position (whether its an overlap or not)
        if i == targetStart:
            mappingOverlapPattern["overlapPattern"].append(chrPotMatchArea[i]["overlap"])
            if chrPotMatchArea[i]["overlap"] == True:
                mappingOverlapPattern["overlap_region_1"]["start"] = i
                mappingOverlapPattern["overlap_region_1"]["end"] = i
                continue

        # Only add overlaps if they are different to the previous overlap
        if chrPotMatchArea[i]["overlap"] != mappingOverlapPattern["overlapPattern"][-1]:
            mappingOverlapPattern["overlapPattern"].append(chrPotMatchArea[i]["overlap"])

        # ensure that the end of the first overlap region is updated
        if chrPotMatchArea[i]["overlap"] == True and i > mappingOverlapPattern["overlap_region_1"]["end"] and mappingOverlapPattern["overlapPattern"] != [True, False, True]:
            mappingOverlapPattern["overlap_region_1"]["end"] = i


        if chrPotMatchArea[i]["overlap"] == True:
            # If this is the first overlap, update the start and end of the first overlap region
            # The first overlap region happens if the False is not in the overlapPattern position list or if the overlap pattern is [False, True]
            if False not in mappingOverlapPattern["overlapPattern"] or mappingOverlapPattern["overlapPattern"] == [False, True]:
                # if the start of the overlap region is 0, update it
                if mappingOverlapPattern["overlap_region_1"]["start"] == 0:
                    mappingOverlapPattern["overlap_region_1"]["start"] = i
                # if the end of the overlap region is less than the current position, update it
                if i > mappingOverlapPattern["overlap_region_1"]["end"]:
                    mappingOverlapPattern["overlap_region_1"]["end"] = i
            # If this is the second overlap, update the start and end of the second overlap region
            # The second overlap region happens if overlap pattern is [True, False, True]
            elif mappingOverlapPattern["overlapPattern"] == [True, False, True]:
                # if the start of the overlap region is 0, update it
                if mappingOverlapPattern["overlap_region_2"]["start"] == 0:
                    mappingOverlapPattern["overlap_region_2"]["start"] = i
                # if the end of the overlap region is less than the current position, update it
                if i > mappingOverlapPattern["overlap_region_2"]["end"]:
                    mappingOverlapPattern["overlap_region_2"]["end"] = i


    printed_info = False
    # Get the updated coordinates for the mapping
    # Don't run loop if the mapping is completely overlapped
    if mappingOverlapPattern["overlapPattern"] == completeOverlapPattern:
        if printed_info == False:
            print("Mapping is completely overlapped. Marked for removal")
            printed_info = True
        return {"queryChr": queryChr, "oldQueryStart": queryStart, "oldQueryEnd": queryEnd, "newQueryStart": newQueryStart, "newQueryEnd": newQueryEnd, "queryStrand": queryStrand, "oldTargetStart": targetStart, "oldTargetEnd": targetEnd, "newTargetStart": newTargetStart, "newTargetEnd": newTargetEnd, "mappingOverlapPattern": mappingOverlapPattern["overlapPattern"], "markForRemoval": True}
    
    # Also don't run loop if the mapping has no overlaps
    if mappingOverlapPattern["overlapPattern"] == noOverlapPattern:
        if printed_info == False:
            print("No overlaps exist for this mapping")
            printed_info = True
        return {"queryChr": queryChr, "oldQueryStart": queryStart, "oldQueryEnd": queryEnd, "newQueryStart": newQueryStart, "newQueryEnd": newQueryEnd, "queryStrand": queryStrand, "oldTargetStart": targetStart, "oldTargetEnd": targetEnd, "newTargetStart": newTargetStart, "newTargetEnd": newTargetEnd, "mappingOverlapPattern": mappingOverlapPattern["overlapPattern"], "markForRemoval": False}

    # Run loop if the mapping has 1 or 2 overlaps
    # If the mapping has 2 overlaps, set the range from the start of the first overlap region to the end of the second overlap region (inclusive)
    if mappingOverlapPattern["overlapPattern"] == [True, False, True]:
        rangeToCheck = list(range(mappingOverlapPattern["overlap_region_1"]["start"], mappingOverlapPattern["overlap_region_1"]["end"] + 1))
        rangeToCheck.extend(list(range(mappingOverlapPattern["overlap_region_2"]["start"], mappingOverlapPattern["overlap_region_2"]["end"] + 1)))
    else:
        rangeToCheck = range(mappingOverlapPattern["overlap_region_1"]["start"], mappingOverlapPattern["overlap_region_1"]["end"] + 1)
    
    # Commence the loop through the overlapping regions
    onSecondOverlapRegion = False
    updatedSecondOverlap = False
    updatedRightOverlap = False
    for i in rangeToCheck:
        # Check if moved from the first overlap region to the second overlap region
        # This will be the case if there is a gap > 1 between the i and the previous i indicating a jump to the next overlap
        if i != rangeToCheck[0] and i - rangeToCheck[rangeToCheck.index(i) - 1] > 1:
            onSecondOverlapRegion = True

        # Check if the current position is an overlap, should be at this point but just an extra precaution
        if not chrPotMatchArea[i]["overlap"] == True:
            sys.exit("Error: The current position is not an overlap region. This should not happen. Please check the code.")
        
        # loop through chrPotMatchArea[i]["chromosomes_matched"] and extract the data from the current chromosome
        for j in chrPotMatchArea[i]["chromosomes_matched"]:
            if j["queryChr"] == queryChr:
                # check which kind of overlap it is
                # Also add +1 to updates because last i will be the last position of the overlap region and you want the first position of the non-overlap region
                if mappingOverlapPattern["overlapPattern"] == twoOverlapPattern:
                    if j["queryMatchedToStrand"] == "-":
                        if printed_info == False:
                            print("Overlap exists on both sides of the mapping, mapping is on the - strand")
                            printed_info = True
                        # Check for which region of overlap you are on
                        # if you are on the first, update the NEW_TARTGET_END and NEW_QUERY_END parameters
                        # if you are on the second, update the NEW_TARTGET_START and NEW_QUERY_START parameters
                        if onSecondOverlapRegion == False:
                            # update the NEW_TARTGET_START and NEW_QUERY_END parameters
                            # Don't add one to the END parameter because the last base of the overlap (where the overlap and the last base of the - strand aligned chromosome is) is 1 base higher than the actaul end of the mapping and the END columns of .paf files are non-inclusive
                            newTargetStart = chrPotMatchArea[i]["target_position"] + 1
                            newQueryEnd = j["locationOnQuery"] 
                        elif onSecondOverlapRegion == True and updatedSecondOverlap == False:
                            # update the NEW_TARTGET_END and NEW_QUERY_START parameters
                            # ONLY update once because the start position will be the 1 base the START of the overlap region (ie. the first loop)
                            # Don't subtract one from the END parameter because the first base of the overlap is one base higher than the actual start of the mapping and the END columns of .paf files are non-inclusive
                            # But add one to the START parameter because the first base of the overlap is one base higher than the actual start of the mapping and the START columns of .paf files are 0 - indexed and inclusive
                            newTargetEnd = chrPotMatchArea[i]["target_position"]
                            newQueryStart = j["locationOnQuery"] + 1
                            updatedSecondOverlap = True
                        
                    elif j["queryMatchedToStrand"] == "+":
                        if printed_info == False:
                            print("Overlap exists on both sides of the mapping, mapping is on the + strand")
                            printed_info = True
                        # Check for which region of overlap you are on
                        # if you are on the first, update the NEW_TARTGET_START and NEW_QUERY_START parameters
                        # if you are on the second, update the NEW_TARTGET_END and NEW_QUERY_END parameters
                        # ONLY update once because the end position will be the 1 base before the START of the overlap region, which in this case is equal to the start of the overlap region because the END columns of .paf files are non-inclusive
                        if onSecondOverlapRegion == False:
                            newTargetStart = chrPotMatchArea[i]["target_position"] + 1 
                            newQueryStart = j["locationOnQuery"] + 1
                        elif onSecondOverlapRegion == True and updatedSecondOverlap == False:
                            newTargetEnd = chrPotMatchArea[i]["target_position"]
                            newQueryEnd = j["locationOnQuery"]
                            updatedSecondOverlap = True
                    
                elif mappingOverlapPattern["overlapPattern"] == leftOverlapPattern:
                    if j["queryMatchedToStrand"] == "-":
                        if printed_info == False:
                            print("Overlap exists on left side, mapping is on - strand")
                            printed_info = True
                        # update the NEW_TARTGET_START and NEW_QUERY_END parameters
                        # Don't add one to the END parameter because the last base of the overlap (where the overlap and the last base of the - strand aligned chromosome is) is 1 base higher than the actaul end of the mapping and the END columns of .paf files are non-inclusive
                        newTargetStart = chrPotMatchArea[i]["target_position"] + 1
                        newQueryEnd = j["locationOnQuery"] 
                        
                    elif j["queryMatchedToStrand"] == "+":
                        if printed_info == False:
                            print("Overlap exists on left side, mapping is on + strand")
                            printed_info = True
                        # update the NEW_TARTGET_START and NEW_QUERY_START parameters
                        newTargetStart = chrPotMatchArea[i]["target_position"] + 1
                        newQueryStart = j["locationOnQuery"] + 1
                    
                elif mappingOverlapPattern["overlapPattern"] == rightOverlapPattern and updatedRightOverlap == False:
                    if j["queryMatchedToStrand"] == "-":
                        if printed_info == False:
                            print("Overlap exists on right side, mapping is on - strand")
                            printed_info = True
                        # update the NEW_TARTGET_END and NEW_QUERY_START parameters
                        # ONLY update once because the start position will be the 1 base the START of the overlap region (ie. the first loop)
                        # Don't subtract one to the END parameter because the first base of the overlap is 1 base higher than the actaul end of the mapping and the END columns of .paf files are non-inclusive
                        newTargetEnd = chrPotMatchArea[i]["target_position"] 
                        newQueryStart = j["locationOnQuery"] + 1
                        updatedRightOverlap = True
                    elif j["queryMatchedToStrand"] == "+":
                        if printed_info == False:
                            print("Overlap exists on right side, mapping is on + strand")
                            printed_info = True
                        # update the NEW_TARTGET_END and NEW_QUERY_END parameters
                        # ONLY update once because the end position will be the 1 base the START of the overlap region (ie. the first loop)
                        # Don't subtract one to the END parameter because the first base of the overlap is 1 base higher than the actaul end of the mapping and the END columns of .paf files are non-inclusive
                        newTargetEnd = chrPotMatchArea[i]["target_position"] 
                        newQueryEnd = j["locationOnQuery"] 
                        updatedRightOverlap = True

                # remove the chromosome from the list
                chrPotMatchArea[i]["chromosomes_matched"].remove(j)
                break

    return {"queryChr": queryChr, "oldQueryStart": queryStart, "oldQueryEnd": queryEnd, "newQueryStart": newQueryStart, "newQueryEnd": newQueryEnd, "queryStrand": queryStrand, "oldTargetStart": targetStart, "oldTargetEnd": targetEnd, "newTargetStart": newTargetStart, "newTargetEnd": newTargetEnd, "mappingOverlapPattern": mappingOverlapPattern["overlapPattern"], "markForRemoval": False}

#  Commence the second loop
print("Starting the second loop to remove overlaps")
finalTable = []

for rowi, row in enumerate(first_pass_table):
    # Print progress
    if rowi % 10 == 0:
        print(f"Progress: {rowi}/{len(first_pass_table)}")

    queryChr = row[0]
    queryStart = int(row[2])
    queryEnd = int(row[3])
    targetStart = int(row[7])
    targetEnd = int(row[8])
    matchingChromosome = row[5].split("_")[-1]
    
    chrPotMatchArea = mappingInfoOnChromosomes[matchingChromosome]
    
    print(f"Checking for overlaps in mapping from chromosome/sequence {queryChr}...")
    new_data = removeOverlap(targetStart, targetEnd, queryChr, queryStart, queryEnd, chrPotMatchArea)

    row[2] = new_data["newQueryStart"]
    row[3] = new_data["newQueryEnd"]
    row[7] = new_data["newTargetStart"]
    row[8] = new_data["newTargetEnd"]

    if new_data["markForRemoval"] == False:
        finalTable.append(row)
        print(f"Removed gaps from a mapping from chromosome/sequence {queryChr} and added to final table") 
    

print("Done removing overlaps from mappings")


# Build the final table string
print("Formatting final table for output...")
outputString = ""
for row in finalTable:
    for i, value in enumerate(row):
        if i!= len(row)-1:
            outputString += str(value) + "\t"
        else:
            outputString += str(value)
    outputString += "\n"
print("Done formatting final table for output")
    
# Output the final table
with open(outFile, 'w') as f:
    print(f"Writing final table to {outFile}...")
    f.write(outputString)
    print("Done")

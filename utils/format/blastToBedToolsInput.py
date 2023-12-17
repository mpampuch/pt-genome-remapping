import argparse
import sys

# Just takes in all the arguments and sets up the help text
parser = argparse.ArgumentParser(description='Updates gff file')
parser.add_argument("-bo", "--blastOutput", help = "Input blast output file. Format should be blast outfmt 6")
parser.add_argument("-of", "--outputFile", help = "Provide a name for an output file (will create this file automatically). This will contain your main output")

args = parser.parse_args()

separator = "\t"

try:
    blastFile = args.blastOutput
    outFile = args.outputFile
except:
    print("Please provide all arguments")
    sys.exit(1)

featureNames = []

# Create dictionary of keys
with open(blastFile,'r') as f:
    lines = f.readlines();
    for line in lines:
        lineVal = line.rstrip()
        parts = lineVal.split(separator)
        featureNames.append(parts[0].rstrip())

dict = {k: [] for k in featureNames}

# Function to add lines from a file to one of the index of the dictionary's value
def addFileLinesToDict(fileName):
    with open(fileName, 'r') as f:
        lines = f.readlines();
        for line in lines:
            lineVal = line.rstrip()
            parts = lineVal.split(separator)
            keyToCheck = parts[0]

            if keyToCheck in dict:
                dict[keyToCheck].append(lineVal)

# Contains our final output
finalOutputString = ""

# Add the blast file lines as values for existing keys in the dictionary
addFileLinesToDict(blastFile)

def breakLine(line, sep):
    lineVal = line.rstrip()
    parts = lineVal.split(sep)
    return parts

for key in dict.keys():
    for index, line in enumerate(dict[key]):
        parts = breakLine(line, separator)
        
        chromosome = str(parts[1])
        start = str(parts[8])
        end = str(parts[9])
        
        # If start > end then swap for bedtools to work properly
        if int(start) > int(end):
            temp = start
            start = end
            end = temp
        
        name = parts[0] + "@" + str(index)
        filler = "0"
        strandInfo = "+"
        
        if parts[12]=="minus":
            strandInfo = "-"
            
        result = chromosome + "\t" + start + "\t" + end + "\t" + name + "\t" + filler + "\t" + strandInfo
        finalOutputString += result + "\n"     
        
# Write to output files
with open(outFile, 'w') as f:
    f.write(finalOutputString)

import argparse
import sys
from Bio import SeqIO
import tempfile
import subprocess
from Bio.Emboss.Applications import NeedleCommandline

# Just takes in all the arguments and sets up the help text
parser = argparse.ArgumentParser(description='Updates gff file')
parser.add_argument("-nf", "--nameFile", help = "Input file containing names of features")
parser.add_argument("-gf", "--geneFastaFile", help = "Input fasta file for feature sequences.")
parser.add_argument("-sf", "--subjectFastaFile", help = "Input fasta file for subject hit sequences. Assumes the fasta headers looks like this (notice the @ placement): >IDgene_Phatr3_J12807_biotypeprotein_coding_descriptionPredictedprotein[Source_UniProtKBTrEMBL3BAcc_B7G006]_gene_idahref3DgenesPhatr3_J12807Phatr3_J12807a_logic_nameensemblgenomes___64478::9:220834-222271(+)@0::9:220835-222271(+)")
parser.add_argument("-bf", "--blastFile", help = "Input the blast file to update. Make sure the LAST column of this file contains the row number.")
parser.add_argument("-obf", "--outBlastFile", help = "Provide the name of the output blast file.")

args = parser.parse_args()

separator = "\t"

try:
    nameFile = args.nameFile
    geneFasta = args.geneFastaFile
    subjectFasta = args.subjectFastaFile
    blastFile = args.blastFile
    outputBlastFile = args.outBlastFile
except:
    print("Please provide all arguments")
    sys.exit(1)

# Input the fasta sequences

def fastaToArray(fastaFile):
    array = []
    fastaObject = SeqIO.parse(open(fastaFile),'fasta')
    for record in fastaObject:
        array.append(record)
    return array    

gene_fasta_sequences = fastaToArray(geneFasta)
subject_fasta_sequences = fastaToArray(subjectFasta)

featureNames = []

# Create array of feature (gene) names
with open(nameFile,'r') as f:
    lines = f.readlines();
    for line in lines:
        name = line.rstrip()
        featureNames.append(name)

geneDict = {k: [] for k in featureNames}
subjectDict = {k: [] for k in featureNames}
 
# Function used for debugging
def printFirstX(yourDict, numberOfItemsToPrint):
    for index, key in enumerate(yourDict.keys()):
        if index > numberOfItemsToPrint:
            break
        print(len(yourDict[key]))
    
    
def buildDict(fastaSeqs, yourDict, splitHeader=False):
    # Build dictonary
    for fasta in fastaSeqs:
        header = fasta.id
        
        if splitHeader==True:
            # Want to remove the extra non-matching (compared to our feature list) part of the header
            parts = header.split("@")
            keyToCheck = parts[0]
        else:
            keyToCheck = header

        if keyToCheck in yourDict:
            yourDict[keyToCheck].append(fasta)
            

buildDict(gene_fasta_sequences, geneDict)      
buildDict(subject_fasta_sequences, subjectDict, True)    
            
#printFirstX(geneDict, 100)
#print("#######################")
#printFirstX(subjectDict,100)
            
def buildFastaOutputString(array):
    result = ""
    # loop over all fasta records in each group
    for index, fastaRecord in enumerate(array):
        # Build a single fasta record
        result += ">"
        result += fastaRecord.id
        result += "\n"
        result += str(fastaRecord.seq)
        result += "\n"
    return result

def countMatchingCharsBetweenStrings(string1, string2):
    numMatching = 0
    for index, char in enumerate(string1):
        if string1[index] is string2[index]:
            numMatching += 1
    return numMatching;

def returnMatchingFirstLastXBases(n, arrayOfStrings):
    firstString = arrayOfStrings[0]
    secondString = arrayOfStrings[1]
    
    firstXCharsOfStringOne = firstString[0:n]
    firstXCharsOfStringTwo = secondString[0:n]
            
    lastXCharsOfStringOne = firstString[-n:]
    lastXCharsOfStringTwo = secondString[-n:]

    numberMatchingFirstXBases = countMatchingCharsBetweenStrings(firstXCharsOfStringOne, firstXCharsOfStringTwo)
    numberMatchingLastXBases = countMatchingCharsBetweenStrings(lastXCharsOfStringOne, lastXCharsOfStringTwo)
    
    result = [numberMatchingFirstXBases, numberMatchingLastXBases]
    return result

def parseNeedleOutput(needleOutput, headerIdentifier, alignSepChar, commentCharacters):
    records = []
    
    lines = needleOutput.splitlines()
    
    recordNum = -1
    oldRecordNum = recordNum
    inHeader = False
    
    header = []
    alignmentSeqs = []
    
    for lineNum, line in enumerate(lines):
        
        cleanedLine = line.rstrip()

        # Only run this on lines that are not whitespace 
        if cleanedLine:
            # Read lines till you hit the header identifier, if we havent already reached it, set inHeader to true
            # If that variable is already true (then we want to exit the header region) so set it to false
            if cleanedLine == headerIdentifier and not inHeader:
                recordNum += 1
                inHeader = True
            elif cleanedLine == headerIdentifier and inHeader:
                inHeader = False

            enteredNewHeader = oldRecordNum != recordNum and recordNum != 0
            doneFinalRecord = lineNum == len(lines) - 1
            
            # Line used to debug the parsing
            #print(str(cleanedLine) + " " + str(oldRecordNum) + " " + str(recordNum) + " " + str(inHeader) + " " + str(enteredNewHeader) + " " + str(doneFinalRecord))
            
            # If we have just entered a new header
            # Then append that finished record to records
            if enteredNewHeader or doneFinalRecord:
                oldRecordNum = recordNum
                newRecord = [header, alignmentSeqs]
                records.append(newRecord)
                header = []
                alignmentSeqs = [] 

            if inHeader and cleanedLine != headerIdentifier:
                    header.append(cleanedLine)
            elif not inHeader and cleanedLine != headerIdentifier:
                firstCharOfLine = cleanedLine[0]
                
                # Will let us seperate the two alignments from each other later. Happens only when we run into ">" and will ignore cases like ">>" or ">>>"
                if firstCharOfLine == alignSepChar and cleanedLine.count(alignSepChar) == 1:
                        alignmentSeqs.append(alignSepChar)

                if firstCharOfLine not in commentCharacters:
                    alignmentSeqs.append(cleanedLine)
                    
    return records

def seperateAlignmentsToArray(alignmentArray, alignSepChar):
    result = []

    alingmentString = ''.join(alignmentArray)
    result = alingmentString.split('>')
    
    # Remove empty parts of the array needed because the string input before splitting was like this >SEQUENCE1>SEQUENCE2
    # Could also just remove the first part of the result array but this is more resistant to any changes we might make
    result = list(filter(None,result))
    
    return result

def getScoreFromHeader(index, header):
    scoreLine = header[index]
    result = scoreLine.replace("# Score: ", "")
    return result

def extractDesiredInfoFromNeedle(recordsArray, alignSepChar):
    result = []
    for record in recordsArray:
        header = record[0]
        alignmentString = record[1]
        
        alignmentInfoParts = seperateAlignmentsToArray(alignmentString, alignSepChar)
        
        matchingFirstAndLastBases = returnMatchingFirstLastXBases(10, alignmentInfoParts)

        score = getScoreFromHeader(12, header)
        newRecord = [score, matchingFirstAndLastBases[0], matchingFirstAndLastBases[1]]
        
        result.append(newRecord)
    return result


alignmentScoresDict = {}

lengthOfSubjectDict = len(subjectDict.keys())

commentCharacters = ["#", ">>>", ">", ">>", ";"]
alignSepChar = ">"
headerIdentifier = "#======================================="

for index, key in enumerate(subjectDict.keys()):
    stage = str(index) + " of " + str(lengthOfSubjectDict) 
    consoleLogString = "Processing " + stage + " : " + key + " through needle alignment"
    print(consoleLogString)
    if len(subjectDict[key]) != 0:
     
        # Make a temporary file
        geneFinalOutputString = buildFastaOutputString(geneDict[key])
        subjectFinalOutputString = buildFastaOutputString(subjectDict[key])
        
        geneTempFile = tempfile.NamedTemporaryFile(suffix='.fasta', prefix='gene')
        subjectTempFile = tempfile.NamedTemporaryFile(suffix='.fasta', prefix='subject')

        with open(geneTempFile.name, 'w') as f:
            f.write(geneFinalOutputString)
        
        with open(subjectTempFile.name, 'w') as f:
            f.write(subjectFinalOutputString)

        command = NeedleCommandline(asequence=geneTempFile.name, bsequence=subjectTempFile.name, gapopen=10, gapextend=0.5, auto = True, aformat = "markx10", stdout=True)
        
        commandString = str(command)
        
        commandParts = commandString.split(" ")
        
        result = subprocess.run(commandParts, stdout=subprocess.PIPE).stdout.decode('utf-8')
 
        # They are stored in the same order as are the records in the subjectDict[key] array of records
        records = parseNeedleOutput(result, headerIdentifier, alignSepChar, commentCharacters)

        # These look like this:
        # records = [ [score index0, first 10 matching index0, last 10 matching index0], [score index1, first 10 matching index1, last 10 matching index1], etc..... ]
        records = extractDesiredInfoFromNeedle(records, alignSepChar)

        for index, record in enumerate(subjectDict[key]):
            header = record.id
            headerParts = header.split("@")
            blastRowNumber = headerParts[1]
            informationForLine = records[index]
            alignmentScoresDict[blastRowNumber] = informationForLine

        geneTempFile.close()
        subjectTempFile.close()


# Console message
print("########Finished global alignments, starting to update blast file##########")


# Now read in the blast file and make a dictonary with the keys being the row number and the values being the entire blast lines
blastDict = {}

with open(blastFile,'r') as f:
    lines = f.readlines();
    for line in lines:
        cleanedLine = line.rstrip()
        partsOfLine = cleanedLine.split("\t")
        rowNum = partsOfLine[len(partsOfLine)-1]
        blastDict[rowNum] = cleanedLine

# Build final output
finalOutPutString = ""
for rowNumber in blastDict.keys():
    rowString = blastDict[rowNumber]
    if rowNumber in alignmentScoresDict:
        # Score
        rowString += "\t"
        rowString += str(alignmentScoresDict[rowNumber][0])
        
        # First 10 matching
        rowString += "\t"
        rowString += str(alignmentScoresDict[rowNumber][1])
        
        # Last 10 matching
        rowString += "\t"
        rowString += str(alignmentScoresDict[rowNumber][2])
        
    rowString += "\n"
    finalOutPutString += rowString
    
with open(outputBlastFile, 'w') as f:
    f.write(finalOutPutString)

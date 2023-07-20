import argparse
import sys


# Just takes in all the arguments and sets up the help text
parser = argparse.ArgumentParser(description='Adds whole chromosome data to the GFF file with remapped features.')
parser.add_argument("-i", "--gffFile", help = "Input GFF file")
parser.add_argument("-o", "--outputFile", help = "Output file")
parser.add_argument("-ns", "--noSave", action = "store_true", help = "Don't save pkl files (default is to save)")
parser.add_argument("-nl", "--noLoad", action = "store_true", help = "Don't load pkl files (default is to load)")

args = parser.parse_args()

separator = "\t"

try:
    inputGff = args.gffFile
    outFile = args.outputFile
except:
    print("Please provide all arguments")
    sys.exit(1)
    
chromosome_info = [["chromosome_1","OU594942","chromosome",1,2608419,".",".",".","ID=chromosome:1;Alias=OU594942.1,chr_1"],
["chromosome_2","OU594942","chromosome",1,2499861,".",".",".","ID=chromosome:2;Alias=OU594943.1,chr_2"],
["chromosome_3","OU594942","chromosome",1,2064744,".",".",".","ID=chromosome:3;Alias=OU594944.1,chr_3"],
["chromosome_4","OU594942","chromosome",1,1629129,".",".",".","ID=chromosome:4;Alias=OU594945.1,chr_4"],
["chromosome_5","OU594942","chromosome",1,1555020,".",".",".","ID=chromosome:5;Alias=OU594946.1,chr_5"],
["chromosome_6","OU594942","chromosome",1,1417157,".",".",".","ID=chromosome:6;Alias=OU594947.1,chr_6"],
["chromosome_7","OU594942","chromosome",1,1124623,".",".",".","ID=chromosome:7;Alias=OU594948.1,chr_7"],
["chromosome_8","OU594942","chromosome",1,1122386,".",".",".","ID=chromosome:8;Alias=OU594949.1,chr_8"],
["chromosome_9","OU594942","chromosome",1,1108211,".",".",".","ID=chromosome:9;Alias=OU594950.1,chr_9"],
["chromosome_10","OU594942","chromosome",1,1107389,".",".",".","ID=chromosome:10;Alias=OU594951.1,chr_10"],
["chromosome_11","OU594942","chromosome",1,1087446,".",".",".","ID=chromosome:11;Alias=OU594952.1,chr_11"],
["chromosome_12","OU594942","chromosome",1,1052234,".",".",".","ID=chromosome:12;Alias=OU594953.1,chr_12"],
["chromosome_13","OU594942","chromosome",1,959323,".",".",".","ID=chromosome:13;Alias=OU594954.1,chr_13"],
["chromosome_14","OU594942","chromosome",1,898576,".",".",".","ID=chromosome:14;Alias=OU594955.1,chr_14"],
["chromosome_15","OU594942","chromosome",1,897230,".",".",".","ID=chromosome:15;Alias=OU594956.1,chr_15"],
["chromosome_16","OU594942","chromosome",1,860830,".",".",".","ID=chromosome:16;Alias=OU594957.1,chr_16"],
["chromosome_17","OU594942","chromosome",1,803256,".",".",".","ID=chromosome:17;Alias=OU594958.1,chr_17"],
["chromosome_18","OU594942","chromosome",1,759571,".",".",".","ID=chromosome:18;Alias=OU594959.1,chr_18"],
["chromosome_19","OU594942","chromosome",1,716929,".",".",".","ID=chromosome:19;Alias=OU594960.1,chr_19"],
["chromosome_20","OU594942","chromosome",1,709265,".",".",".","ID=chromosome:20;Alias=OU594961.1,chr_20"],
["chromosome_21","OU594942","chromosome",1,629758,".",".",".","ID=chromosome:21;Alias=OU594962.1,chr_21"],
["chromosome_22","OU594942","chromosome",1,587839,".",".",".","ID=chromosome:22;Alias=OU594963.1,chr_22"],
["chromosome_23","OU594942","chromosome",1,557589,".",".",".","ID=chromosome:23;Alias=OU594964.1,chr_23"],
["chromosome_24","OU594942","chromosome",1,546843,".",".",".","ID=chromosome:24;Alias=OU594965.1,chr_24"],
["chromosome_25","OU594942","chromosome",1,516884,".",".",".","ID=chromosome:25;Alias=OU594966.1,chr_25"],
["mitochondrion","Genbank","region",1,77356,".","+",".","ID=HQ840789.1:1..77356;Dbxref=taxon:2850;Is_circular=true;Name=MT;gbkey=Src;genome=mitochondrion;mol_type=genomic DNA"],
["chloroplast","ASM15095v2","chromosome",1,117369,".",".",".","ID=chromosome:chloroplast;Alias=EF067920.1,chr_chloroplast"]]


def buildGffTableWithChromosomeInfo(gffFile, chromosome_info, separator):
    chr_info_table = chromosome_info
    gffTable = []
    with open(gffFile,'r') as f:
        lines = f.readlines();
        currentChr = ""
        
        for i, line in enumerate(lines):
          lineVal = line.rstrip()
          parts = lineVal.split(separator)
          
          # If this is the first line, add the chromosome info to the new gff table
          if i == 0 or not currentChr == parts[0] :
            gffTable.append(chr_info_table[0])
            # update the current chromosome
            currentChr = chr_info_table[0][0]
            # remove the appended chromosome info from the chr_info_table
            chr_info_table.pop(0)
          
          # Add the line to the new gff table
          # Create if condition to prevent double adding mitochondrion and chloroplast information
          if parts[8] not in ["ID=chromosome:chloroplast;Alias=EF067920.1,chr_chloroplast", "IDchromosome_chloroplast_AliasEF067920.1chr_chloroplast___69248", "ID=HQ840789.1:1..77356;Dbxref=taxon:2850;Is_circular=true;Name=MT;gbkey=Src;genome=mitochondrion;mol_type=genomic DNA", "IDHQ840789.1_1..77356_Dbxreftaxon_2850_Is_circulartrue_NameMT_gbkeySrc_genomemitochondrion_mol_typegenomicDNA___69354"]:
            gffTable.append(parts)
            
          # print progress
          if i % 10000 == 0:
            print(f"Added line {i} out of {len(lines)} to new gff table")
            
    return gffTable
  
final_gff_table_list = buildGffTableWithChromosomeInfo(inputGff, chromosome_info, separator)
print("Done building new gff table")
print("Configuring table for output")
outputString = ""
for i, row in enumerate(final_gff_table_list):
    for j, value in enumerate(row):
        if j!= len(row)-1:
            outputString += str(value) + "\t"
        else:
            outputString += str(value)
    outputString += "\n"
    
    # Print progress
    if i % 10000 == 0:
        print(f"Configured line {i} out of {len(final_gff_table_list)} for output")

with open(outFile, 'w') as f:
    print(f"Writing final table to {outFile}...")
    f.write(outputString)
    print("Done")
    
"""
All the data that this script was created with

# data obtained from the utils/chr_sizes/get_sizes.sh script
# New assembly GenBank accession ID: OU594942
chromosome_1	OU594942	chromosome	1	2608419	.	.	.	ID=chromosome:1;Alias=OU594942.1,chr_1
chromosome_2	OU594942	chromosome	1	2499861	.	.	.	ID=chromosome:2;Alias=OU594943.1,chr_2
chromosome_3	OU594942	chromosome	1	2064744	.	.	.	ID=chromosome:3;Alias=OU594944.1,chr_3
chromosome_4	OU594942	chromosome	1	1629129	.	.	.	ID=chromosome:4;Alias=OU594945.1,chr_4
chromosome_5	OU594942	chromosome	1	1555020	.	.	.	ID=chromosome:5;Alias=OU594946.1,chr_5
chromosome_6	OU594942	chromosome	1	1417157	.	.	.	ID=chromosome:6;Alias=OU594947.1,chr_6
chromosome_7	OU594942	chromosome	1	1124623	.	.	.	ID=chromosome:7;Alias=OU594948.1,chr_7
chromosome_8	OU594942	chromosome	1	1122386	.	.	.	ID=chromosome:8;Alias=OU594949.1,chr_8
chromosome_9	OU594942	chromosome	1	1108211	.	.	.	ID=chromosome:9;Alias=OU594950.1,chr_9
chromosome_10	OU594942	chromosome	1	1107389	.	.	.	ID=chromosome:10;Alias=OU594951.1,chr_10
chromosome_11	OU594942	chromosome	1	1087446	.	.	.	ID=chromosome:11;Alias=OU594952.1,chr_11
chromosome_12	OU594942	chromosome	1	1052234	.	.	.	ID=chromosome:12;Alias=OU594953.1,chr_12
chromosome_13	OU594942	chromosome	1	959323	.	.	.	ID=chromosome:13;Alias=OU594954.1,chr_13
chromosome_14	OU594942	chromosome	1	898576	.	.	.	ID=chromosome:14;Alias=OU594955.1,chr_14
chromosome_15	OU594942	chromosome	1	897230	.	.	.	ID=chromosome:15;Alias=OU594956.1,chr_15
chromosome_16	OU594942	chromosome	1	860830	.	.	.	ID=chromosome:16;Alias=OU594957.1,chr_16
chromosome_17	OU594942	chromosome	1	803256	.	.	.	ID=chromosome:17;Alias=OU594958.1,chr_17
chromosome_18	OU594942	chromosome	1	759571	.	.	.	ID=chromosome:18;Alias=OU594959.1,chr_18
chromosome_19	OU594942	chromosome	1	716929	.	.	.	ID=chromosome:19;Alias=OU594960.1,chr_19
chromosome_20	OU594942	chromosome	1	709265	.	.	.	ID=chromosome:20;Alias=OU594961.1,chr_20
chromosome_21	OU594942	chromosome	1	629758	.	.	.	ID=chromosome:21;Alias=OU594962.1,chr_21
chromosome_22	OU594942	chromosome	1	587839	.	.	.	ID=chromosome:22;Alias=OU594963.1,chr_22
chromosome_23	OU594942	chromosome	1	557589	.	.	.	ID=chromosome:23;Alias=OU594964.1,chr_23
chromosome_24	OU594942	chromosome	1	546843	.	.	.	ID=chromosome:24;Alias=OU594965.1,chr_24
chromosome_25	OU594942	chromosome	1	516884	.	.	.	ID=chromosome:25;Alias=OU594966.1,chr_25
# Organelle information with Original Metadata
chloroplast	ASM15095v2	chromosome	1	117369	.	.	.	ID=chromosome:chloroplast;Alias=EF067920.1,chr_chloroplast
mitochondrion	Genbank	region	1	77356	.	+	.	ID=HQ840789.1:1..77356;Dbxref=taxon:2850;Is_circular=true;Name=MT;gbkey=Src;genome=mitochondrion;mol_type=genomic DNA
# Organelle information with Cleaned up Metadata
chloroplast	ASM15095v2	chromosome	1	117369	.	.	.	IDchromosome_chloroplast_AliasEF067920.1chr_chloroplast___69248
mitochondrion	Genbank	region	1	77356	.	+	.	IDHQ840789.1_1..77356_Dbxreftaxon_2850_Is_circulartrue_NameMT_gbkeySrc_genomemitochondrion_mol_typegenomicDNA___69354	```

# New chromosome Aliases

# Molecule name	GenBank sequence
# Chromosome 1	OU594942.1
# Chromosome 2	OU594943.1
# Chromosome 3	OU594944.1
# Chromosome 4	OU594945.1
# Chromosome 5	OU594946.1
# Chromosome 6	OU594947.1
# Chromosome 7	OU594948.1
# Chromosome 8	OU594949.1
# Chromosome 9	OU594950.1
# Chromosome 10	OU594951.1
# Chromosome 11	OU594952.1
# Chromosome 12	OU594953.1
# Chromosome 13	OU594954.1
# Chromosome 14	OU594955.1
# Chromosome 15	OU594956.1
# Chromosome 16	OU594957.1
# Chromosome 17	OU594958.1
# Chromosome 18	OU594959.1
# Chromosome 19	OU594960.1
# Chromosome 20	OU594961.1
# Chromosome 21	OU594962.1
# Chromosome 22	OU594963.1
# Chromosome 23	OU594964.1
# Chromosome 24	OU594965.1
# Chromosome 25	OU594966.1

# Information with Original Metadata
1	ASM15095v2	chromosome	1	2535400	.	.	.	ID=chromosome:1;Alias=CM000605.1,chr_1
2	ASM15095v2	chromosome	1	1497954	.	.	.	ID=chromosome:2;Alias=CM000606.1,chr_2
3	ASM15095v2	chromosome	1	1460046	.	.	.	ID=chromosome:3;Alias=CP001142.1,chr_3
4	ASM15095v2	chromosome	1	1360148	.	.	.	ID=chromosome:4;Alias=CM000607.1,chr_4
5	ASM15095v2	chromosome	1	1098047	.	.	.	ID=chromosome:5;Alias=CM000608.1,chr_5
6	ASM15095v2	chromosome	1	1035082	.	.	.	ID=chromosome:6;Alias=CM000609.1,chr_6
7	ASM15095v2	chromosome	1	1029019	.	.	.	ID=chromosome:7;Alias=CM000610.1,chr_7
8	ASM15095v2	chromosome	1	1007773	.	.	.	ID=chromosome:8;Alias=CM000611.1,chr_8
9	ASM15095v2	chromosome	1	1002813	.	.	.	ID=chromosome:9;Alias=CM000612.1,chr_9
10	ASM15095v2	chromosome	1	976485	.	.	.	ID=chromosome:10;Alias=CM000613.1,chr_10
11	ASM15095v2	chromosome	1	945026	.	.	.	ID=chromosome:11;Alias=CP001141.1,chr_11
12	ASM15095v2	chromosome	1	901853	.	.	.	ID=chromosome:12;Alias=CM000614.1,chr_12
13	ASM15095v2	chromosome	1	887524	.	.	.	ID=chromosome:13;Alias=CM000615.1,chr_13
14	ASM15095v2	chromosome	1	829358	.	.	.	ID=chromosome:14;Alias=CM000616.1,chr_14
15	ASM15095v2	chromosome	1	814910	.	.	.	ID=chromosome:15;Alias=ABQD01000058.1,chr_15
16	ASM15095v2	chromosome	1	764225	.	.	.	ID=chromosome:16;Alias=CM000618.1,chr_16
17	ASM15095v2	chromosome	1	703943	.	.	.	ID=chromosome:17;Alias=ABQD01000062.1,chr_17
18	ASM15095v2	chromosome	1	702471	.	.	.	ID=chromosome:18;Alias=ABQD01000063.1,chr_18
19	ASM15095v2	chromosome	1	690427	.	.	.	ID=chromosome:19;Alias=CM000621.1,chr_19
20	ASM15095v2	chromosome	1	683011	.	.	.	ID=chromosome:20;Alias=CM000622.1,chr_20
21	ASM15095v2	chromosome	1	662217	.	.	.	ID=chromosome:21;Alias=CM000623.1,chr_21
22	ASM15095v2	chromosome	1	591336	.	.	.	ID=chromosome:22;Alias=CM000624.1,chr_22
23	ASM15095v2	chromosome	1	512847	.	.	.	ID=chromosome:23;Alias=CM000625.1,chr_23
24	ASM15095v2	chromosome	1	511739	.	.	.	ID=chromosome:24;Alias=CM000626.1,chr_24
25	ASM15095v2	chromosome	1	497271	.	.	.	ID=chromosome:25;Alias=CM000627.1,chr_25
26	ASM15095v2	chromosome	1	441226	.	.	.	ID=chromosome:26;Alias=CM000628.1,chr_26
27	ASM15095v2	chromosome	1	404298	.	.	.	ID=chromosome:27;Alias=CM000629.1,chr_27
28	ASM15095v2	chromosome	1	387582	.	.	.	ID=chromosome:28;Alias=ABQD01000085.1,chr_28
29	ASM15095v2	chromosome	1	384256	.	.	.	ID=chromosome:29;Alias=CM000631.1,chr_29
30	ASM15095v2	chromosome	1	317207	.	.	.	ID=chromosome:30;Alias=CM000632.1,chr_30
31	ASM15095v2	chromosome	1	258242	.	.	.	ID=chromosome:31;Alias=ABQD01000095.1,chr_31
32	ASM15095v2	chromosome	1	157053	.	.	.	ID=chromosome:32;Alias=CM000634.1,chr_32
33	ASM15095v2	chromosome	1	87967	.	.	.	ID=chromosome:33;Alias=CM000635.1,chr_33
chloroplast	ASM15095v2	chromosome	1	117369	.	.	.	ID=chromosome:chloroplast;Alias=EF067920.1,chr_chloroplast
mitochondrion	Genbank	region	1	77356	.	+	.	ID=HQ840789.1:1..77356;Dbxref=taxon:2850;Is_circular=true;Name=MT;gbkey=Src;genome=mitochondrion;mol_type=genomic DNA

# Information with Cleaned up Metadata
1	ASM15095v2	chromosome	1	2535400	.	.	.	IDchromosome_1_AliasCM000605.1chr_1___0
2	ASM15095v2	chromosome	1	1497954	.	.	.	IDchromosome_2_AliasCM000606.1chr_2___27908
3	ASM15095v2	chromosome	1	1460046	.	.	.	IDchromosome_3_AliasCP001142.1chr_3___44000
4	ASM15095v2	chromosome	1	1360148	.	.	.	IDchromosome_4_AliasCM000607.1chr_4___49587
5	ASM15095v2	chromosome	1	1098047	.	.	.	IDchromosome_5_AliasCM000608.1chr_5___53104
6	ASM15095v2	chromosome	1	1035082	.	.	.	IDchromosome_6_AliasCM000609.1chr_6___55859
7	ASM15095v2	chromosome	1	1029019	.	.	.	IDchromosome_7_AliasCM000610.1chr_7___58647
8	ASM15095v2	chromosome	1	1007773	.	.	.	IDchromosome_8_AliasCM000611.1chr_8___61404
9	ASM15095v2	chromosome	1	1002813	.	.	.	IDchromosome_9_AliasCM000612.1chr_9___63904
10	ASM15095v2	chromosome	1	976485	.	.	.	IDchromosome_10_AliasCM000613.1chr_10___6831
11	ASM15095v2	chromosome	1	945026	.	.	.	IDchromosome_11_AliasCP001141.1chr_11___9266
12	ASM15095v2	chromosome	1	901853	.	.	.	IDchromosome_12_AliasCM000614.1chr_12___11843
13	ASM15095v2	chromosome	1	887524	.	.	.	IDchromosome_13_AliasCM000615.1chr_13___14200
14	ASM15095v2	chromosome	1	829358	.	.	.	IDchromosome_14_AliasCM000616.1chr_14___16660
15	ASM15095v2	chromosome	1	814910	.	.	.	IDchromosome_15_AliasABQD01000058.1chr_15___18674
16	ASM15095v2	chromosome	1	764225	.	.	.	IDchromosome_16_AliasCM000618.1chr_16___20887
17	ASM15095v2	chromosome	1	703943	.	.	.	IDchromosome_17_AliasABQD01000062.1chr_17___22797
18	ASM15095v2	chromosome	1	702471	.	.	.	IDchromosome_18_AliasABQD01000063.1chr_18___24508
19	ASM15095v2	chromosome	1	690427	.	.	.	IDchromosome_19_AliasCM000621.1chr_19___26262
20	ASM15095v2	chromosome	1	683011	.	.	.	IDchromosome_20_AliasCM000622.1chr_20___31821
21	ASM15095v2	chromosome	1	662217	.	.	.	IDchromosome_21_AliasCM000623.1chr_21___33561
22	ASM15095v2	chromosome	1	591336	.	.	.	IDchromosome_22_AliasCM000624.1chr_22___35133
23	ASM15095v2	chromosome	1	512847	.	.	.	IDchromosome_23_AliasCM000625.1chr_23___36786
24	ASM15095v2	chromosome	1	511739	.	.	.	IDchromosome_24_AliasCM000626.1chr_24___38036
25	ASM15095v2	chromosome	1	497271	.	.	.	IDchromosome_25_AliasCM000627.1chr_25___39100
26	ASM15095v2	chromosome	1	441226	.	.	.	IDchromosome_26_AliasCM000628.1chr_26___40266
27	ASM15095v2	chromosome	1	404298	.	.	.	IDchromosome_27_AliasCM000629.1chr_27___41221
28	ASM15095v2	chromosome	1	387582	.	.	.	IDchromosome_28_AliasABQD01000085.1chr_28___42150
29	ASM15095v2	chromosome	1	384256	.	.	.	IDchromosome_29_AliasCM000631.1chr_29___43222
30	ASM15095v2	chromosome	1	317207	.	.	.	IDchromosome_30_AliasCM000632.1chr_30___47870
31	ASM15095v2	chromosome	1	258242	.	.	.	IDchromosome_31_AliasABQD01000095.1chr_31___48574
32	ASM15095v2	chromosome	1	157053	.	.	.	IDchromosome_32_AliasCM000634.1chr_32___49148
33	ASM15095v2	chromosome	1	87967	.	.	.	IDchromosome_33_AliasCM000635.1chr_33___49417
chloroplast	ASM15095v2	chromosome	1	117369	.	.	.	IDchromosome_chloroplast_AliasEF067920.1chr_chloroplast___69248
mitochondrion	Genbank	region	1	77356	.	+	.	IDHQ840789.1_1..77356_Dbxreftaxon_2850_Is_circulartrue_NameMT_gbkeySrc_genomemitochondrion_mol_typegenomicDNA___69354	```
"""
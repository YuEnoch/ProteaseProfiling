#Yu Enoch 30-07-2024
#Purpose: Takes 5-amino acid mutagenesis sequencing data to derive a merged cleaved peptide counts table, and frequencies of amino acids and nucleotides.

#Please set up parameters.txt and store all FastQ or Gzip files in the same folder as this file
#Parameters: Name of Experiment, Treatments, Control, Number of Treatments, Replicates, Mutagenesis Positions, Row Sum Filter, Gzip, Single/Parallel Reads, File Names

#Specifications:
# 1. Uses single-end reads or assembled files of paired-end reads (e.g., by PEAR)
# 2. Can work on any n-mer, but requires a file with all possible n-mers (e.g., NNK5_combinations.txt)
# 3. Mutagenesis positions are inclusive (51-65 for 15 nucleotides, corresponds to line[50:65])
# 4. The number of replicates are same for experiments/control
# 5. Quality Threshold: removes any read where 2 consecutive bases have a PHRED score below 28

import sys
import os
import subprocess
from subprocess import call


def getDetails(str, type):
    str = str.strip().split(' = ')
    str = str[1]
    if type == "s":
        return str
    elif type == "i":
        return int(str.replace(" ", ""))
    elif type == "l":
        str = str.split(', ')
        return str

print("Starting Pipeline for Sequencing Data: \n\nThis file would filter, trim, and translate the sequence reads of cleaved phages with randomized 5-mers.\nSubsequent figures and processing would be done in R, as denoted in the R markdown files.\nComplete supplementary data of sequences and counts can be located in SRA and GEO databases.")
inp = int(input("\nWould the processing take place for Proteases (1) or Mixtures (2): "))

if inp == 1:
    parameterFile = "parameters_protease.txt" #either parameters_protease.txt or parameters_mixture.txt
elif inp == 2:
    parameterFile = "parameters_mixture.txt"
else:
    parameterFile = "parameters.txt"

print("Pipeline in Progress...")

input = open(parameterFile,'r')

input.readline()
experimentName = getDetails(input.readline(), 's')
experiments = getDetails(input.readline(), 'l')
control = getDetails(input.readline(), 'l')

input.readline()
treatments = getDetails(input.readline(), 'i')
n_plicates = getDetails(input.readline(), 'i')  #duplicates, triplicates, etc.
treatmentNames = getDetails(input.readline(), 'l')
allTreatments = getDetails(input.readline(), 'l')

input.readline()
input.readline()
firstPosition = getDetails(input.readline(), 'i') - 1 #-1 because Python counts 0 as the 1st Position
lastPosition = getDetails(input.readline(), 'i')
n_mer = getDetails(input.readline(), 'i')
if ((lastPosition - firstPosition)/3 != n_mer):
    sys.exit("Error: Mismatch on Positions and N-Mers")
elif (len(treatmentNames) * n_plicates != treatments or len(experiments) + len(control) != len(treatmentNames) or treatments != len(allTreatments)):
    sys.exit("Error: Mismatch on Number of Treatments")
else:
    firstPosition = str(firstPosition)
    lastPosition = str(lastPosition)
    n_mer = str(n_mer)
minCount = str(getDetails(input.readline(), 'i'))

input.readline()
input.readline()
gzipChoice = getDetails(input.readline(), 's') #Are the FastQ files in Gzip? 0 = No, 1 = Yes
if gzipChoice.lower() == 'yes':
    gzipChoice = str(1)
else:
    gzipChoice = str(0)
singleParallel = getDetails(input.readline(), 's').lower()
filenames = getDetails(input.readline(), 'l')

def find_files(fileName):              #Find File Location (if FastQ is in same folder)
    dir_path = os.getcwd()            
    for root, dirs, files in os.walk(dir_path):
        if fileName in files: 
            return root+'/'+str(fileName)    

#clean_fastq.py: takes in raw FastQ data (single-end read) and filters sequencing reads by Quality and Matched Seeds, to obtain good Sense sequences (and Junk data)
#getPeptides.py: takes in Sense sequences and translates the Mutagenesis Region to obtain a table of Peptide Frequencies

for i in range(len(allTreatments)):
    name = allTreatments[i]
    fileName1 = find_files(filenames[i])
    if singleParallel == 'single':
        fileName2 = 'SingleReadOnly'
    else:
        fileName2 = find_files(fileNameFormat(i, 2, gzipChoice))
    print(name, fileName1, fileName2)

    call(["python", "clean_fastq.py", name, fileName1, fileName2, gzipChoice, firstPosition, lastPosition])
    print(name, "clean_fastq.py done")
    call(["python", "getPeptides.py", name, firstPosition, lastPosition])
    print(name, "getPeptides.py done")

    call(["python", "getAminoAcidCount.py", name, n_mer])
    call(["python", "getNucleicAcidCount.py", name, firstPosition, lastPosition])
    call(["python", "ObservedExpectedPlot.py", name, n_mer])

#Generates Observed/Expected, Nucleotide, Amino Acid Frequencies for one treatment (combining its n-plicates)
for i in range(int(len(allTreatments)/n_plicates)): 
    treatmentList = ""
    for j in range(i*n_plicates, i*n_plicates+n_plicates):
        treatmentList+= allTreatments[j] + ','
    treatmentList = treatmentList[:-1]
    treatment = treatmentNames[i]
    print(treatmentList)
    call(["python", "MergeCounts.py", treatment, str(n_plicates), treatmentList, str(0), n_mer])
    call(["python", "MergeCounts.py", treatment, str(n_plicates), treatmentList, str(1), n_mer])
    call(["python", "ObservedExpectedPlot.py", treatment, n_mer])

#Generates a merged Obs/Exp Table for all treatments
call(["python", "FrequencyAnalysis.py", ','.join(treatmentNames), ','.join(allTreatments)])

print("Merging Peptide Data")

#Merge_for_DESeq2.py: merges Peptide Frequencies across all treatments into one file, for subsequent Deseq2 analysis
treatmentList = ','.join(allTreatments)
call(["python", "Merge_for_DESeq2.py", experimentName, str(treatments), treatmentList, minCount, n_mer])

#Generates a complete counts table of cleaved peptides
#Subsequent Analysis is done via R scripts (DeepProteaseProfiling.Rmd and PreliminaryAnalysis.Rmd)

print("Pipeline Complete")



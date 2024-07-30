#YuEnoch 30-07-2024
#Merge_for_Deseq2.py 

#Purpose: takes in all Peptide Counts for each treatment and merges all into one Table (ready to be analyzed by DESeq2), with peptides of low counts filtered out
#Changes across Experiments: please alter on parameters.txt
# 1. Location of Mutagenesis
# 2. Mutagenesis for 5 Amino Acids
# 3. Number of Treatments

import sys
import re
import mod575
import gzip
import os
import csv
from collections import defaultdict

#Opens the necessary files for input and output
experimentName = sys.argv[1]
numberOfTreatments = sys.argv[2] #Example: 15
num = int(numberOfTreatments)
treatments = sys.argv[3] #Example: C1, C2, C3, E1, E2, E3, H1, H2, H3, V1, V2, V3, N1, N2, N3
treatments = treatments.split(',')
minCount = int(sys.argv[4]) #Minimum Count for Each Peptide (across the row), rowsum filter
n_mer = int(sys.argv[5])

combinations=open('NNK' + str(n_mer) + '_combinations.txt','r')	#All possible combinations for n-mer amino acids
outfile=open(experimentName+"_data.txt",'w')
outfileCSV=open(experimentName+"_data.csv",'w', newline='')
writer = csv.writer(outfileCSV)
header = ["AA_seq"] + treatments
writer.writerow(header)

noreads=0
peptides=dict() #Empty Peptides Dictionary
for line in combinations:	
    line=line.strip()
    if line in peptides:	
        print('non_unique')
    else:				
        peptides[line]=0	


for i in range(num):            #Creates a Peptides Dictionary for each treatment, then goes through every Peptides File and adds to that Dictionary
    input1 = open(treatments[i]+'_peptide', 'r')
    globals()['peptide%s' % i] = peptides.copy()
    for line in input1:			
        row=line.strip().split() 	#separates into peptide, peptide count, sequence distribution
        peptide=row[0] 			
        peptide_count=row[1]
        if peptide in globals()['peptide%s' % i]:
            globals()['peptide%s' % i][peptide]=peptide_count
        else:
            print('sth is wrong with dictionary!')        
    input1.close()
        
for peptide in peptides:
    sum = 0
    for i in range(num):
        sum += int(globals()['peptide%s' % i][peptide])
    if sum<minCount: 	#Threshold. If total counts for that peptide is below minCount, it is not added to final result
        noreads+=1
        continue
    else:
        list = [peptide]
        for i in range(num):
            list.append(str(globals()['peptide%s' % i][peptide]))
        print ('\t'.join(list), file=outfile)
        writer.writerow(list)
outfile.close()
outfileCSV.close()


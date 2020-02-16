# This is modified from https://github.com/levisimons/oysterGut/blob/master/demultiplex.py
#This script will take two input files (fwd and rev reads), trim primer sequences, and demultiplex into per sample reads.

import os,re,sys
import gzip

lineNum=1

print("Starting...")

#Create barcode and index libraries. 
barcode515F = ['TCAGC','GTATC','GCTAC','ACGCA','CGCGT','CTGGT','GCGTT','GGAAC','AAGCC','AGCTT','CCCTT','CGCAC','GGTGT','GGCAG','TGATA','TGTGC']
# barcode is found in line 2 of fwd fastq (_1.fq) file (characters 5-9)
index926R = ['CGTGAT', 'ACATCG','GCCTAA','TGGTCA','CACTGT','ATTGGC','GATCTG','TCAAGT','TGACAT','GGACGG','GCGGAC','TTTCAC','CCGGTG','ATCGTG','TGAGTG','CGCCTG']
# reverse complement of reverse index
#index926R = ['TCACG','GATGT','TAGGC','GACCA','CAGTG','CCAAT','AGATC','CTTGA','TGTCA','CGTCC','TCCGC','TGAAA','ACCGG', 'ACGAT','ACTCA','AGGCG']
# index is found in line 1 of either fastq file (last 6 characters)

print("Index libraries loaded.")

fwdFileSuffix = '_1.fastq'
revFileSuffix = '_2.fastq'
UndeterminedFwdFileName = 'Undetermined',str(fwdFileSuffix)
UndeterminedFwdFileName = ''.join(UndeterminedFwdFileName)
UndeterminedRevFileName = 'Undetermined',str(revFileSuffix)
UndeterminedRevFileName = ''.join(UndeterminedRevFileName)
fwdIndex = 'Blank'
revIndex = 'Blank'

#Read in fwd and rev reads
fwdRead = open('small_MO_1_paired.fq')
revRead = open('small_MO_2_paired.fq')

print("Fwd and Rev read files loaded.")
print("Attempting to demultiplex...")

for fwdLine, revLine in zip(fwdRead, revRead):
    if lineNum%4 == 1:
        #Read in the sequence identifier.  This the first of four lines. Extract reverse index.
        fwdIdentifier=str(fwdLine)
        revIdentifier=str(revLine)
        revIndex=revLine[-9:-3]
        #revIndex = revIndex.replace('\\n','')
        #revIndex = revIndex.rstrip("\n\r")
        #print(fwdIdentifier)
        #print(revIdentifier)
        #print(revIndex)
    if lineNum%4 == 2:
        # Read in the sequence
        fwdIndex=fwdLine[4:9]
        fwdSequence=fwdLine[28:]
        #print(fwdLine)
        #print(fwdSequence)
        #print(fwdIndex)
        revSequence=revLine[20:]
    if lineNum%4 == 3:
        fwdQI = str(fwdLine)
        revQI = str(revLine)
    if lineNum%4 == 0:
        #Read in quality scores
        fwdQuality = fwdLine[28:]
        revQuality = revLine[20:]
        #make sure that fwd and rev indicies are in the index library. if yes, make a new file
        if fwdIndex in barcode515F:
            if revIndex in index926R:
                fwdIndexNum = barcode515F.index(fwdIndex)+1
                revIndexNum = index926R.index(revIndex)+1
                fwdFileName = 'F',str(fwdIndexNum),'R',str(revIndexNum),str(fwdFileSuffix)
                fwdFileName = ''.join(fwdFileName)
                revFileName = 'F',str(fwdIndexNum),'R',str(revIndexNum),str(revFileSuffix)
                revFileName = ''.join(revFileName)
                fwdOutput = open(fwdFileName,'a')
                fwdOutputLine = fwdIdentifier, fwdSequence, fwdQI, fwdQuality
                fwdOutputLine = ''.join(fwdOutputLine)
                fwdOutput.write(fwdOutputLine)
                revOutput = open(revFileName,'a')
                revOutputLine = revIdentifier, revSequence, revQI, revQuality
                revOutputLine = ''.join(revOutputLine)
                revOutput.write(revOutputLine)
                fwdIndex = 'Blank'
                revIndex = 'Blank'
            #if rev index not in library, add read to undetermined file
            else:
                UndeterminedFwdOutput = open(UndeterminedFwdFileName,'a')
                UndeterminedFwdOutputLine = fwdIdentifier, fwdSequence, fwdQI, fwdQuality
                UndeterminedFwdOutputLine = ''.join(UndeterminedFwdOutputLine)
                UndeterminedFwdOutput.write(UndeterminedFwdOutputLine)
                UndeterminedRevOutput = open(UndeterminedRevFileName,'a')
                UndeterminedRevOutputLine = revIdentifier, revSequence, revQI, revQuality
                UndeterminedRevOutputLine = ''.join(UndeterminedRevOutputLine)
                UndeterminedRevOutput.write(UndeterminedRevOutputLine)
                fwdIndex = 'Blank'
                revIndex = 'Blank'
        #if fwd index not in library, add read to undetermined file
        else:
            UndeterminedFwdOutput = open(UndeterminedFwdFileName,'a')
            UndeterminedFwdOutputLine = fwdIdentifier, fwdSequence, fwdQI, fwdQuality
            UndeterminedFwdOutputLine = ''.join(UndeterminedFwdOutputLine)
            UndeterminedFwdOutput.write(UndeterminedFwdOutputLine)
            UndeterminedRevOutput = open(UndeterminedRevFileName,'a')
            UndeterminedRevOutputLine = revIdentifier, revSequence, revQI, revQuality
            UndeterminedRevOutputLine = ''.join(UndeterminedRevOutputLine)
            UndeterminedRevOutput.write(UndeterminedRevOutputLine)
            fwdIndex = 'Blank'
            revIndex = 'Blank'
    #print(lineNum)
    lineNum=lineNum+1
        

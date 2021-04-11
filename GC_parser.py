###################################################################
#Description: this program imports the Fasta file: GCF_005508785.1_pea_aphid_22Mar2018_4r6ur_genomic.fna (pea aphid genome)
#               The user specifies which chromosome to analyze based on which position of the records used.

#Output: a plot of the GC content across regions (10,000 bases) of the genome
###################################################################

import pandas as pd #use to make dataframes
import matplotlib.pyplot as plt #use to graph GC contents
from Bio import SeqIO
records = list(SeqIO.parse('GCF_005508785.1_pea_aphid_22Mar2018_4r6ur_genomic.fna','fasta'))

import xlrd

import time
start = time.asctime()
print(start)

########################################################################################################################
def GC_calc(sequence):
    #function is passed a string variable (sequence)
    ## the GC content of the sequence will be calculated and returned
    ### store counts of G and C
    GC_count = 0
    for base in sequence:
        if base=="G" or base=="C":
            GC_count+=1
    #calculate proportion that is GC in seq
    GC_prop = GC_count/len(sequence)
    return GC_prop
########################################################################################################################
########################################################################################################################
def slidingWindow(sequence,window_length):
    #function is passed a sequence of entire chromosome
    ## will examine the GC content of a specified sliding window and compare to overall GC content of chromosome
    startPos = 0
    for i in range(round(len(sequence) / window_length)):

        yield sequence[startPos:startPos + window_length]

        startPos = startPos + window_length
########################################################################################################################
########################################################################################################################
########################################################################################################################
#   change Chromosome here
## position 0: chromosomeX, position 1: chromosomeA1, position 2: chromosomeA2, position 3: chromosomeA3
########################################################################################################################
########################################################################################################################
########################################################################################################################
sequence = records[2]
print(sequence.description)
sequence = str(sequence.seq) #store the sequence as a string
#remove 'N' characters & making all character uppercase
sequence = sequence.upper().replace('N','')


print(str(len(sequence)))

#define the overall GC content of the current chromosome
GC_content = GC_calc(sequence)
print('past overall GC content')
#loop through every sub-sequence in the sliding window range
## make a list storing the GC conent of the sub-region & starting position

currentStartPos=0 #stores the current starting position on in the loop below
subSeq_GC = []
########################
for sub_seq in slidingWindow(sequence,10000):
    #determine the GC content of the current sub-seq
    sub_GC = GC_calc(sub_seq)
    subSeqEntry = [currentStartPos,sub_GC]


    subSeq_GC.append(subSeqEntry) #add current entry to subSeq_GC list containing found GC contents
    currentStartPos+=1


########################
#create a data frame to graph
df = pd.DataFrame(subSeq_GC,columns=['start position (10,000 base intervals)','GC content'])
########################
########################################################################################################################
########################################################################################################################
# change the title of the plot here
########################################################################################################################
########################################################################################################################
#graph dataframe to see changes in GC content at positions
df.plot(x='start position (10,000 base intervals)',y='GC content',title='chromosomeA2')
plt.show()





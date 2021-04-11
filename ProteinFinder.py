###################################################################
#Description: this program imports the Fasta file: GCF_005508785.1_pea_aphid_22Mar2018_4r6ur_genomic.fna (pea aphid genome)
#               The user specifies which chromosome to analyze based on which position of the records used.
##              The program will run through all posibble reading frames of a chromosome sequence and will split into
###              possible protein sequences in between Stop Codons, ...Seq1...(stopcodon)...Seq2...(stopcodon)...

#Output: A text file it exported with sample descriptions of a protein and number in the order retrieved, followed by a
##          candidate possible protein sequence found

#typical run time is 11 min for an entire chromosome
###################################################################

import time
from Bio import SeqIO
from Bio import Seq
#######################################
#import genome data, grab sequence for chormosome of choice
records = list(SeqIO.parse('GCF_005508785.1_pea_aphid_22Mar2018_4r6ur_genomic.fna','fasta'))
########################################################################################################################
#   change Chromosome here
## position 0: chromosomeX, position 1: chromosomeA1, position 2: chromosomeA2, position 3: chromosomeA3
########################################################################################################################
chromosome = records[1]
sequence = chromosome.seq
print('description'+str(chromosome.description))
print('total size of chromosome: '+str(len(sequence)))
#######################################
start = time.asctime()
print(start)

########################################################################################################################
def checkStopCodon(seq):
    #give a sequence at will return if it ENDs with a stop codon
    stopCodons = ['TAA', 'TAG', 'TGA']

    for SC in stopCodons: #loop through stopcodons
        if seq.endswith(SC):
            return True #if stopcodon return True
        else:
            pass
    return False #if no stop codon was matched

########################################################################################################################
########################################################################################################################
# sequence = 'ATGTTAGTTGATAG'
candidateSequences = []
current_seq=''
counter=0 #keep track of where at in algorithm

#loop through seq for all possible reading frame and retrieve sequence betweeen stop codons that are larger than 60 AAs (180 base pairs)
for i in range(3):
    #### check if current_seq from the end of last reading frame had a stop codon before resetting
    SC_status = checkStopCodon(current_seq)  # give function current codon to check if stopcodon
    if SC_status:
        if len(current_seq) > 180:  # and determine size to see if store seq IF is at least 60 AAs (180 Bases)
            # translate seq to AA sequence

            current_seq = Seq.translate(current_seq)
            candidateSequences.append(current_seq)
    ###################
    #loop through each reading frame possible
    current_seq = ''
    current_codon = ''
    print('NEXT reading frame')

    for pos in range(len(sequence)-i): #loop through every base in the range of the sequence (minus current reading frame)
        ##########################
        counter+=1
        if counter%100000000==0:
            current_time=time.asctime()
            print('100,000,000 bases checked')
            print(current_time)
        ##########################
        ##########################

        #check if completed a codon

        if len(current_codon)<3: # if not complete codon add another

            current_codon+= str(sequence[pos+i]).upper() #add i to pos to keep track of right current reading frame

            if len(current_codon)==3:
                # print('current seq: '+current_seq)
                # print('current codon: '+ current_codon)
            #if codon IS complete check if is a stop codon
            #    loop through all stopCodons

                SC_status = checkStopCodon(current_codon) #give function current codon to check if stopcodon
                #returns True if found stop codon

                if SC_status: #if find a stop codon at the end of the current seq then stop
                    # print('stop codon found: ' +current_codon)
                    current_seq+=current_codon #add codon to current seq and check on length

                    if len(current_seq)>180: #and determine size to see if store seq IF is at least 60 AAs (180 Bases)
                        # translate seq to AA sequence
                        current_seq = Seq.translate(current_seq)
                        candidateSequences.append(current_seq)
                        current_seq=''#reset seq
                        current_codon=''#start the current codon at the current pos in i before it moves to next


                    else:
                        #if is too small then don't store but reset current_seq & move on down the sequence
                        current_seq=''
                        current_codon=''#start the current codon at the current pos in i before it moves to next


                else:
                #no stop codon found
                    # print('adding codon to seq '+current_codon +'  '+current_seq)
                    current_seq+=current_codon #add to current seq
                    current_codon='' #start the current codon at ''
###################################################################################
#### check if current_seq from the end of the last reading frame had a stop codon before resetting
SC_status = checkStopCodon(current_seq)  # give function current codon to check if stopcodon
if SC_status:
    candidateSequences.append(current_seq)
###################
print('number of candidate sequences: '+str(len(candidateSequences)))
###################
count=1
with open('Candidate_proteins.txt','w') as f:
    for i in candidateSequences:
        header='> protein '+str(count)
        f.write("%s\n" %header)
        f.write("%s\n" %i)
        count+=1
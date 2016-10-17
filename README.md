'''
Andras Palfi
9/27/2016


A python implementation of the overlap alignment problem used in 
bioinformatics sequence alignment.

'''
import sys
from collections import defaultdict
import numpy
from Bio import SeqIO

#Import the FAFSA file from memory. I used an external library called Biopython to do this so
#I could worry more about the overlap allignment problem rather than manual string cleanup
def importFASTA(file):
	
	handle = open(file, "rU" )
	records=list(SeqIO.parse(handle, "fasta"))
	handle.close()

	input_seq_1=records[0].seq
	input_seq_2=records[-1].seq

	return input_seq_1, input_seq_2

'''

Build the matrix used for the overlap problem 
this function takes as input two sequences of 
DNA strings, which were obtained in the previous
method using Biopython. 

I first converted the 'Seq' object into a string,
and we also save its length, which is useful later.

I made two copies of the same matrix one for the path and other for
score
We return all of these elements to be used in the 
overlap allign method


'''
def build_matrix(sequence1,sequence2):

	seq1string=str(sequence1)
	seq2string=str(sequence2)

	seq1Length=len(seq1string)
	seq2Length=len(seq2string)

	MATRIX = [[0 for j1 in	 xrange(seq2Length+1)] for i1 in xrange(seq1Length+1)]
	MATRIX_2 = [[0 for j2 in xrange(seq2Length+1)] for i2 in xrange(seq1Length+1)]
	
	return MATRIX, MATRIX_2, seq1Length, seq2Length, seq1string, seq2string





def overlap_allign(fillMatrix,checkMatrix, dna_length_1, dna_length_2, sequence_string_1, sequence_string_2):
	

#The double for loop used to score the data. We use xrange for speed. Using normal range()
# is extremley slow with this large dataset as it literally builds out all of the values 
#whereas xrange does lazy evaluation. THis is still slow, but a bit faster then normal range
	for i in xrange(1, dna_length_1+1):
		for j in xrange(1, dna_length_2+1):
			
			minusBOTH=fillMatrix[i-1][j-1]
			minusI=fillMatrix[i-1][j]
			minusJ=fillMatrix[i][j-1]		
			
			score1= minusBOTH + [1,-2][sequence_string_1[i-1] == sequence_string_2[j-1]]
			score2= minusI - 2
			score3= minusJ +1
			

			score=[]
			score.append(score1)
			score.append(score2)
			score.append(score3)
			checkMatrix[i][j] = score.index(max(score))
			
	seq1ALLIGN = sequence_string_1[:i] 
	seq2ALLIGN =  sequence_string_2[:j]

	if i == dna_length_1:
		iMAX=i 

	if j == dna_length_2:
		jMAX=j

	score_max=fillMatrix[iMAX][jMAX]


#Here is where I shift the strings and try to replace indels with underscores. 
	while i*j != 0 :
		if checkMatrix[i][j] == 1:
			i = i - 1
			seq2ALLIGN = seq2ALLIGN[:j] + '_' + seq2ALLIGN[:j]
		elif checkMatrix[i][j] == 2:
			j = j - 1
			seq1ALLIGN = seq1ALLIGN[:i] + '_' + seq1ALLIGN[:i]
		else:
			i = i -1
			j  = j -1

	
	return str(score_max), seq1ALLIGN, seq2ALLIGN

#The main method were I call all of the functions from and print the results in 
#the desired format of:
#Score
#Allignment 1
#Allignment 2
def main():
	seq1,seq2=importFASTA('example.fasta')
	m1,m2, sl1, sl2, ss1, ss2=build_matrix(seq1,seq2)
	results = overlap_allign(m1,m2, sl1, sl2, ss1, ss2)

	for result in results:
		print(result)

if __name__ == '__main__':
    sys.exit(main())

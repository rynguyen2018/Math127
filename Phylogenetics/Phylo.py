import random
import numpy as np 
import math
import argparse
import sys

def JCdistance(seq1, seq2):
	differences=0
	if len(seq1)!= len(seq2):
		print("HEY SCREW YOU! I DON'T KNOW WHERE YOU LEARNED YOUR PHYLOGENETICS FROM, BUT, HERE IN AMERICA, WE USE SAME LENGTH SEQUENCES FOR PHYLOGENETICS")
		sys.exit()
	else:
		for base1, base2 in zip(seq1,seq2):
			if base1!= base2:
				differences+=1 
	#return differences/len(seq1)
	return -0.75*math.log(1- ((4/3)*differences/len(seq1)))


def main(): 	
	parser = argparse.ArgumentParser(description='Process some sequence data and make a phylogenetic tree out of it.')
	parser.add_argument('-i','--input', help= "reads in file with sequences in FASTA format")

	args= parser.parse_args()
	file= args.input
	
	num_seq=0
	sequence_array=[]
	with open(file) as seqfile:
		for line in seqfile:
			if ">" in line: 
				num_seq+=1
			else:
				sequence_array.append(line.rstrip())

	JC_array= np.empty([num_seq,num_seq])
	for i in range(0, num_seq):
		for j in range(0, num_seq):
			JC_array[i][j]= JCdistance(sequence_array[i], sequence_array[j])
	print(JC_array)
	
	JC_array_perturb= JC_array
	for _ in range(0,4):
		x=random.randint(0,num_seq-1)
		y=random.randint(0,num_seq-1) 
		if x!=y:
			JC_array_perturb[x][y]=np.add(JC_array[x][y], np.random.normal(0,0.2))

	np.save("JC_array", JC_array)
	np.save("JC_array_perturb", JC_array_perturb)

if __name__ == '__main__':
	main()










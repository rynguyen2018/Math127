import random
import numpy as np 
import math
import argparse
import sys




Markov_array= np.array([[0.25,0.25,0.25,0.25 ],[0.25,0.25,0.25,0.25 ],[0.25,0.25,0.25,0.25 ],[0.25,0.25,0.25,0.25 ]])

def mutate(base_A, base_T, base_C, base_G, original_base): 
	A_mutate= base_A
	T_mutate= base_A+base_T
	C_mutate= base_A+base_T+ base_G

	rand_num= random.random()
	if rand_num<=base_A:
		return "A"
	elif rand_num>base_A and rand_num<=T_mutate:
		return "T"
	elif rand_num>T_mutate and rand_num<= C_mutate:
		return "C"
	return "G"

def JCdistance(seq1, seq2):
	differences=0
	if len(seq1)!= len(seq2):
		print("HEY SCREW YOU! I DON'T KNOW WHERE YOU LEARNED YOUR PHYLOGENETICS FROM, BUT, HERE IN AMERICA, WE USE SAME LENGTH SEQUENCES FOR PHYLOGENETICS")
		sys.exit()
	else:
		for base1, base2 in zip(seq1,seq2):
			if base1!= base2:
				differences+=1 
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

	np.save("JC_array", JC_array)

if __name__ == '__main__':
	main()

# point_mutation=random.randint(0, len(seq1)-1)
# print(point_mutation)
# #for base in seq1:
# seq2= seq1[0:point_mutation] +mutate(0.25,0.25,0.25,0.25, seq1[point_mutation]) +seq1[point_mutation+1: len(seq1)+1]

# print(seq1, seq2)
# print(JCdistance(seq1,seq2))










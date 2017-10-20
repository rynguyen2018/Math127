import random
import numpy as np 
import math
import argparse
import sys
import re


def mutate(Markov_array, original_base):
	mutate_dict= {'A':0, 'T':1, 'C': 2, 'G':3} 
	mutate_x_coordinate=mutate_dict[original_base] 

	A_mutate= Markov_array[mutate_x_coordinate][0]
	T_mutate= A_mutate+Markov_array[mutate_x_coordinate][1]
	C_mutate= T_mutate+Markov_array[mutate_x_coordinate][2]

	rand_num= random.random()
	if rand_num<=A_mutate:
		return "A"
	elif rand_num>A_mutate and rand_num<=T_mutate:
		return "T"
	elif rand_num>T_mutate and rand_num<= C_mutate:
		return "C"
	return "G"

def createTree(seq1,time):
	if time==0:
		return seq1
	else:
		point_mutation1=random.randint(0, len(seq1)-1)
		point_mutation2=random.randint(0, len(seq1)-1)

		seq2= seq1[0:point_mutation1] + mutate(Markov_array, seq1[point_mutation1]) + seq1[point_mutation1+1: len(seq1)+1]
		seq3= seq1[0:point_mutation2] + mutate(Markov_array, seq1[point_mutation2]) + seq1[point_mutation2+1: len(seq1)+1]

		kid2= createTree(seq2,time-1)
		kid3= createTree(seq3, time-1)
		return '(' + kid2+ ',' + kid3 + ')'
alpha=0.8

Markov_array= np.array([[1-alpha,alpha/3,alpha/3,alpha/3 ],
						[alpha/3,1-alpha,alpha/3,alpha/3 ],
						[alpha/3,alpha/3,1-alpha,alpha/3 ],
						[alpha/3,alpha/3,alpha/3,1- alpha ]])

seq1= 'ATGCTGCTCGCTCATGAATTG'


tree= createTree(seq1, 3)
sequence_list= tree.replace("(","")
sequence_list= sequence_list.replace(")","")
sequence_list= sequence_list.split(",")
perfect_tree= open("perfect_data_seq.fasta", "w")
for i in range(1, len(sequence_list)+1):
	title=">"+"seq"+str(i)
	perfect_tree.write( title+ "\n")
	perfect_tree.write(sequence_list[i-1]+ "\n")
perfect_tree.close()

for i in range(0, len(sequence_list)):
	tree= re.sub(sequence_list[i], "seq" +str(i), tree)

tree_data= open("perfect_data.nwk","w")
tree_data.write(tree+';')
tree_data.close()







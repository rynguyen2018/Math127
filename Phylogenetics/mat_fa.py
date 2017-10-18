from scipy.io import loadmat
###turns data into fasta format because I'm not a savage and I strictly follow fasta format
hiv_data = loadmat("flhivdata.mat")

fasta_file= open("HIV_data.fasta", "w")
min_seq_length=99999999999
title_list=[]
seq_list=[] 
for key, value in hiv_data.items(): 
	if "_" not in key: 
		title= ">"+key +"\n"
		title_list.append(title)
		#fasta_file.write(title)
		seq= value[0].replace("'","")
		if len(seq)< min_seq_length:
			min_seq_length= len(seq)
		seq_list.append(seq)
fasta_file= open("HIV_data.fasta", "w")
for i in range(0, len(title_list)): 
	fasta_file.write(title_list[i])
	sequence= seq_list[i][:min_seq_length]+ "\n"
	fasta_file.write(sequence)

fasta_file.close()
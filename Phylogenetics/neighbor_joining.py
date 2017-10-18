import numpy as np 
import sys
#array of dissimilarity scores, delta
delta_array= np.array([[0,0.31,1.01,0.75,1.03],
             [0.31,0,1,0.69,0.90],
             [1.01,1,0,0.61,0.42],
             [0.75,0.69,0.61,0,0.37],
             [1.03,0.9,0.42,0.37,0]])

num_species= len(delta_array[0])

def getM_vals(delta_array, num_species):
	R_values= []

	"""the R values are the sum of all dissimilarity scores pertaining to a particular species i. 
	This is simply the sum of all values in rows""" 
	for i in range(0, num_species): 
		R_values.append(sum(delta_array[i]))
	Mij_array= np.empty([num_species, num_species])
	for row in range(0, num_species):
		for col in range(0, num_species):
			if row==col: 
				Mij_array[row][col]=0 
			else: 
				Mij_array[row][col]= (num_species-2)*delta_array[row][col] - R_values[row]- R_values[col]
	return Mij_array

def min_M(Mij_array, count_array):
	min_row, min_col= -1, -1
	min_value= 999999999999999
	# find minimum value of Mij_array 
	for x in range(0, len(Mij_array[0])):
		for y in range(0, len(Mij_array[:,0])):
			if Mij_array[x][y]< min_value:
				min_row=x
				min_col=y
				min_value= Mij_array[x][y]
	return min_row, min_col

def make_tree(Mij_array, delta_array, num_species, count_array): 
	if num_species==3:
		###End game

		
		#print(ThreePointFormula(delta_array,1,2))
	else: 	
		node1, node2= min_M(Mij_array, count_array)
		new_delta_array, distance1, distance2= ThreePointFormula(delta_array, node1, node2, num_species)
		print(count_array[node1]+1,": " ,distance1, count_array[node2]+1,": ", distance2)
		count_array= [count_array[taxa] for taxa in range(0, len(count_array)) if taxa!= node1 and taxa!= node2]

		new_dimensions= [_ for _ in range(0, len(new_delta_array))]
		new_dimensions.remove(node1)
		new_dimensions.remove(node2)
		ixgrid= np.ix_(new_dimensions, new_dimensions)
		new_delta_array= new_delta_array[ixgrid]
		new_Mij_array= getM_vals(new_delta_array, num_species-1)
		make_tree(new_Mij_array, new_delta_array, num_species-1, count_array)


def ThreePointFormula(delta_array, node1, node2, num_species):
	distanceiv, distancejv= getClusterDistance(node1, node2, delta_array, num_species)

	new_delta_array= np.zeros([len(delta_array[0])+1, len(delta_array[0])+1])
	new_delta_array[:len(delta_array[0]), :len(delta_array[0])]= delta_array
	
	for row in range(0, len(new_delta_array[0])):
		for col in range(0, len(new_delta_array[0])):
			if row==col: 
				new_delta_array[row][col]= 0
			elif col== len(delta_array[0]) and row!= node1 and row != node2 and col!= node1 and col!= node2:
				new_delta_array[row][col]= 0.5*(delta_array[node1][row]+ delta_array[node2][row]- (distanceiv+distancejv))
			elif row== len(delta_array[0]):
				new_delta_array[row][col]= new_delta_array[col][row]
	return new_delta_array, distanceiv, distancejv

def getClusterDistance(node1, node2, delta_array, num_species): 
	distance_node1=0
	distance_node2=0
	for x in range(0, len(delta_array[0])):
		if x != node1 and x!=node2: 
			distance_node1+= delta_array[node1][x]
			distance_node2+= delta_array[node2][x]

	distance_node1, distance_node2= distance_node1/(num_species-2), distance_node2/(num_species-2)
	distance1, distance2= 0.5*(delta_array[node1][node2] + distance_node1-distance_node2), 0.5*(delta_array[node1][node2] + distance_node2-distance_node1)
	return distance1, distance2

Mij_array= getM_vals(delta_array, num_species)
newdistance= make_tree(Mij_array, delta_array, num_species, [_ for _ in range(0, num_species)])
print(newdistance)
from matplotlib import pyplot as plt
from project_fst import project_fst
import numpy as np
import mrcfile
import sys
import scipy.ndimage

a = np.array([1,0,0])
b = np.array([0,1,0])
c = np.cross(a,b)

R1 = np.transpose(np.array([a,b,c]))
R2 = np.transpose(np.array([b,a,np.cross(b,a)]))

zika_file = mrcfile.open('zika_153.mrc')
rho = zika_file.data
	
image1 = project_fst(rho, R1)
image2 = project_fst(rho, R2)

N = image1.shape[0]

def getLine(image, theta, num_points = 100):
	x = np.arange(image.shape[1])
	y = np.arange(image.shape[0])

	if theta == 0:  
		x0 = 0
		y0 = (N - 1)/2  
		x1 = N 
		y1 = (N - 1)/2  
	elif (theta <= np.pi/4 and theta >= -np.pi/4 ) or (theta >= 3*np.pi/4 and theta <= 5*np.pi/4 ):
	 	x0 = 0
	 	y0 = np.tan(theta) * (N - 1)/2 + ((N - 1)/2)
	 	x1 = N
	 	y1 = (N-1)/2 - np.tan(theta) * (N-1)/2
	else:
	 	x0 = np.tan(theta - np.pi/2 ) * (N - 1)/2 + ((N - 1)/2) 
	 	y0 = N
	 	x1 = (N-1)/2 - np.tan(theta - np.pi/2) * (N-1)/2
	 	y1 = 0

	x_values = np.linspace(x0, x1, num_points)
	y_values = np.linspace(y0, y1, num_points)

	z_values = scipy.ndimage.map_coordinates(image, np.vstack((x_values,y_values)))
	#fig, axes = plt.subplots(nrows=2)
	#axes[0].imshow(image)
	#axes[0].plot([x0, x1], [y0, y1], 'ro-')
	#axes[0].axis('image')
	#axes[1].plot(z_values)
	#plt.show()
	return z_values

def commonLine(image1, image2): 
	commonLine_array = []
	commonLine_product_array = []
	for theta_i in range(360):
		zi = getLine(image1, (theta_i * np.pi *2)/360)
		for theta_j in range(360):
			zj = getLine(image2, (theta_j * np.pi *2)/360)
			commonLine_array.append(zi)
			commonLine_array.append(zj)
			commonLine_product_array.append(np.dot(zi,zj))
	max_index, max_value = max(enumerate(commonLine_product_array), key=lambda p: p[1])
	return max_index, max_value, commonLine_array

def getArcLengths(cos_theta):
	return np.arccos(cos_theta/4)

def getGlobalAngles(A,B,C): 
	return np.arccos((np.cos(C)- np.cos(A)*np.cos(N))/(np.sin(A) * np.sin(B))

h, i, j = commonLine(image1,image2)

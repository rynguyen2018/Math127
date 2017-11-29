import numpy as np 
import mrcfile
from scipy.interpolate import RegularGridInterpolator as RGI
from skimage.color import rgb2gray as rg
import sys
from matplotlib import pyplot as plt
import Rotation_Matrix as RM

def project_fst(mol, R):
	""" 
		mol is an NXNXN array that contains the electron density of the protein thing
		R is the rotation matrix that symbolizes the viewing direction of the microscope
	"""

	print( "Now running 3D fourier Transform")

	###perform 3D transform on molecule 
	N = len(mol)
	rho_hat = np.fft.fftn(mol)
	rho_hat = np.fft.fftshift(rho_hat)
	wx, wy, wz = np.meshgrid(np.arange(-(N-1)/2, (N-1)/2 +1, dtype= int),np.arange(-(N-1)/2, (N-1)/2 +1,dtype= int),np.arange(-(N-1)/2, (N-1)/2 +1,dtype= int))
	rho_hat=((-1)**np.abs(wx+wy+wz)/(N**3)*rho_hat) #np.exp(np.pi*1j*(wx+wy+wz))/(N**3))*(rho_hat)
	
	#print("Now running Linear Interpolation")
	N_range = np.arange(-(N-1)/2, (N-1)/2 +1)
	#print("N_range found")
	rho_hat = RGI((N_range,N_range,N_range), rho_hat, bounds_error= False, fill_value=0)

	### make 2D mesh that takes slice of rho hat which NxNx3 grid 
	#print("Creating 2D Sampling grid")
	eta_x, eta_y= np.meshgrid(np.arange(-(N-1)/2, (N-1)/2 +1, dtype= int),np.arange(-(N-1)/2, (N-1)/2 +1, dtype= int), indexing= "ij")

	eta_x= eta_x[...,np.newaxis]
	eta_y= eta_y[...,np.newaxis]
	
	R= np.transpose(R)
	a= R[0]
	b= R[1]

	grid = eta_x*a + eta_y*b

	sample_grid = rho_hat(grid)
	
	em_slice=(-1)**np.abs(eta_x[...,0]+eta_y[...,0])*(N**2)*sample_grid

	em_slice = np.fft.ifftshift(em_slice)
	em_slice = np.fft.ifftn(em_slice)
	#print(np.linalg.norm(np.imag(em_slice)))

	return np.real(em_slice)


#zika_file=mrcfile.open('zika_153.mrc')
#rho = zika_file.data

#print("Now running project_fst")
#R= RM.get_rotation_matrix()
#image = project_fst(rho, R)

#print("Now showing image")
#plt.imshow(image, cmap= "Greys")
#plt.show()
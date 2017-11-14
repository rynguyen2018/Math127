import numpy as np 
import mrcfile
from scipy.interpolate import RegularGridInterpolator as RGI
import sys
from matplotlib import pyplot as plt

def project_fst(mol, R=0):
	""" 
		mol is an NXNXN array that contains samples of the molecule rho
		R is the rotation matrix that symbolizes the viewing direction of the microscope
	"""


	print( "Now running 3D fourier Transform")




	###perform 3D transform on molecule 
	N = mol.shape[0]
	rho_hat = np.fft.fftshift(mol)
	rho_hat = np.fft.fftn(rho_hat)
	print(np.real(rho_hat))
	
	print("Now running Linear Interpolation")
	N_range = np.linspace(-N, N, N)
	rho_hat = RGI((N_range,N_range,N_range), rho_hat, bounds_error= False, fill_value=0)
	#sys.exit()
	### make 2D mesh that takes slice of rho hat which NxNx3 grid 
	### to do this, we must find the etas
	print("Creating 2D Sampling grid")
	eta_x = np.zeros((N,N))
	eta_y = np.zeros((N,N))
	for i in range(N):
		for j in range(N):
			eta_x[i,j]= -N/2 +1 + j
			eta_y[i,j]= -N/2 + 1 +i
	eta_x= eta_x[...,np.newaxis]
	eta_y= eta_y[...,np.newaxis]

	a = np.array([1,0,0])
	b = np.array([0,1,0])
	grid= eta_x*a + eta_y*b
	#print(grid)

	sample_grid= rho_hat(grid)

	em_slice= np.fft.fftshift(sample_grid)
	em_slice= np.fft.ifftn(em_slice)

	return np.real(em_slice)

	### we now sample rho_hat on the grid= rho_hat_f(grid)
	### take the inverse fourier transform of this grid 
	### get the real portion 
	### plt.imshow(image)

	### 

	#mol_size = mol.shape[0]
	#N = np.linspace(0,mol_size-1,mol_size)
	
	#ThreeD_sampling_grid = RegularGridInterpolator((N, N, N), mol, method='linear', bounds_error= False, fill_value=0)


	### Okay so we have done the linear interpolation for the electron densities for this particular molecule mol.
	### However, this still doesn't explain how we get

	#a_vector= [1,0,0]
	#b_vector= [0,1,0]
	#basis = np.array([a_vector, b_vector])
	#TwoD_sampling_grid=np.zeros(N,N)
	#for i in range(0,mol_size):
	#	for j in range(0,mol_size):
	#		for k in range(0, mol_size):
	#			TwoD_sampling_grid[i,k]= np.linalg.solve(basis, ThreeD_sampling_grid([i,j,k]))
	#print(ThreeD_sampling_grid([1,2,3]))
	#sys.exit()
	#eta_array=np.array([[1,0,1],[0,1,2]])



zika_file=mrcfile.open('zika_153.mrc')
rho = zika_file.data
#print(rho)
print("Now running project_fst")
image = project_fst(rho)

print("Now showing image")
plt.imshow(image)
plt.show()
#print(rho.shape[0])
# 	#	2D_sampling_grid= 


	
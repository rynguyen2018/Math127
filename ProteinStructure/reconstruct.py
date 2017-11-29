import numpy as np
from project_fst import project_fst
import mrcfile
import Rotation_Matrix as RM
import numpy as np 
import mrcfile
from scipy.interpolate import RegularGridInterpolator as RGI
import sys


def reconstruct(image_array, R_array): 
	N = image_array[0].shape[0] 
	B = np.zeros((N,N,N)) 
 	
	#print(R_array[0])

	for i in range(0, len(image_array)):
		
		print("Computing image_hat")
		image_hat = np.fft.fftn(image_array[i])
		image_hat = np.fft.fftshift(image_hat)
		freq_range = np.arange(-(N-1)/2, (N-1)/2 +1, dtype= int) 
		wx_j, wy_j = np.meshgrid(freq_range, freq_range)
		image_hat = ((-1)**np.abs(wx_j+wy_j)/(N**2)*image_hat)
		image_hat = np.tile(image_hat[...,np.newaxis], [1,1,N])

		print("Computing l_j_hat")
		#wlx_j, wly_j, wlz_j = np.meshgrid(freq_range, freq_range, freq_range)
		wlz_j = freq_range#np.meshgrid(freq_range)
		l_j = N * np.sinc(N * np.pi * wlz_j)
		l_j_hat = np.tile(l_j[np.newaxis,np.newaxis,...],[N,N,1])
		#l_j_hat = ((-1)**np.abs(wlx_j+wly_j+wlz_j)/(N**3)*l_j_hat)
		
		print("Creating rotated grid")

		#np.array((freq_range,freq_range,3))
		rot_grid = np.meshgrid(freq_range, freq_range, freq_range)
		#d = np.zeros((freq_range, freq_range, freq_range, 3))
		rot_grid = np.tensordot(rot_grid, R_array[i], axes = )
		print(type(rot_grid))

		sys.exit()
		print("Computing b_j_hat")
		
		b_j_hat = image_hat * l_j_hat
		b_j_hat_fnc = RGI((freq_range,freq_range,freq_range), b_j_hat, bounds_error= False, fill_value=0)
		b_j_hat = b_j_hat_fnc(rot_grid)

		b_j_hat = np.fft.ifftshift(b_j_hat)
		b_j = np.fft.ifftn(b_j_hat)

		print(np.real(b_j))
		#b_j = ((-1)**np.abs(wx_j+wy_j+wz_j)/(N**3)*b_j_hat)
		B += np.real(b_j)
	return np.float32(np.real(B))

num_images = 1
zika_file = mrcfile.open('zika_153.mrc')
rho = zika_file.data

image_array = []
R_array = []
for _ in range(0,num_images):
	R= RM.get_rotation_matrix()
	
	image = project_fst(rho, R)
	image_array.append(image)
	R_array.append(R)
print("Now getting image")
back_img = reconstruct(image_array, R_array)
print("Now saving")
with mrcfile.new('full2.mrc', overwrite=True) as mrc:
	mrc.set_data(back_img)

#print(back_img.shape)
#with file('back_img.mrc', 'w') as outfile:
#    for slice_2d in image:
#        np.savetxt(outfile, slice_2d)
#np.savetxt('back_img.mrc', back_img, delimiter=',')
#back_obj = open("back_img.mrc", "w")
#back_obj.write(back_img)
#back_obj.close()


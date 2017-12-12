import numpy as np
from project_fst import project_fst
import mrcfile
import Rotation_Matrix as RM
import numpy as np 
import mrcfile
from scipy.interpolate import RegularGridInterpolator as RGI
import sys
import time

start = time.time()

def reconstruct(image_array, R_array, filter = False): 
	N = image_array[0].shape[0] 
	B = np.zeros((N,N,N)) 
	L = np.zeros((N,N,N))

	for i in range(0, len(image_array)):
		
		print("Computing image_hat")
		image_hat = np.fft.fftn(image_array[i])
		image_hat = np.fft.fftshift(image_hat)
		freq_range = np.arange(-(N-1)/2, (N-1)/2 +1, dtype= int) 
		wx_j, wy_j = np.meshgrid(freq_range, freq_range)
		image_hat = (((-1)**np.abs(wx_j+wy_j))/N**2)*image_hat
		image_hat = np.tile(image_hat[...,np.newaxis], [1,1,N])

		print("Computing l_j_hat")
		wz_j = freq_range
		l_j = N * np.sinc(N * np.pi * wz_j)
		l_j_hat = np.tile(l_j[np.newaxis,np.newaxis,...],[N,N,1])


		print("Creating rotated grid")
		wlx, wly, wlz = np.meshgrid(freq_range, freq_range, freq_range)
		rot_grid = np.zeros((N,N,N,3))

		for j in range(N):
			for k in range(N):
				for l in range(N):
					point = np.array([wlx[j,k,l], wly[j,k,l], wlz[j,k,l]]) 
					rot_grid[j,k,l] = np.dot(R_array[i], point) 
		
		image_hat_fnc = RGI((freq_range,freq_range,freq_range), image_hat, bounds_error= False, fill_value=0)
		l_j_hat_fnc = RGI((freq_range,freq_range,freq_range), l_j_hat, bounds_error= False, fill_value=0)

		image_hat = image_hat_fnc(rot_grid)
		l_j_hat = l_j_hat_fnc(rot_grid)

		print("Now computing B_j")
		b_j_hat = image_hat * l_j_hat
		b_j_hat = ((-1)**np.abs(wlx+wly+wlz)*(N**3))*b_j_hat
		if filter: 
			b_j_hat = np.divide(np.real(b_j_hat), np.real(l_j_hat), out=np.zeros_like(np.real(b_j_hat)), where=np.real(l_j_hat)!=0)
			b_j_hat = np.fft.ifftshift(b_j_hat)
			b_j = np.fft.ifftn(b_j_hat)
			#L += np.real(l_j_hat)
			#B += np.real(b_j_hat)
		else:
			b_j_hat = np.fft.ifftshift(b_j_hat)
			b_j = np.fft.ifftn(b_j_hat)
		B += np.real(b_j)
	#if filter:
		#B = np.divide(B, L, out=np.zeros_like(B), where=L!=0)
		#B = np.fft.ifftshift(B)
		#B = np.fft.ifftn(B)
	return np.float32(np.real(B))
	#return np.float32(np.real(B)/np.real(L))

num_images = 15
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
back_img = reconstruct(image_array, R_array, filter= False)
#back_img[np.isinf(back_img)] = np.float32(0)
#back_img[np.isneginf(back_img)] = np.float32(0)
#back_img[np.isnan(back_img)] = np.float32(0)
#print(back_img)
print("Now saving")
with mrcfile.new('full2.mrc', overwrite=True) as mrc:
	mrc.set_data(back_img)
end = time.time()
print(end - start)
#print(back_img.shape)
#with file('back_img.mrc', 'w') as outfile:
#    for slice_2d in image:
#        np.savetxt(outfile, slice_2d)
#np.savetxt('back_img.mrc', back_img, delimiter=',')
#back_obj = open("back_img.mrc", "w")
#back_obj.write(back_img)
#back_obj.close()


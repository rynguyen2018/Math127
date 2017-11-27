import numpy as np
from project_fst import project_fst
import mrcfile
import Rotation_Matrix as RM
import numpy as np 
import mrcfile
from scipy.interpolate import RegularGridInterpolator as RGI
import sys


def reconstruct(image, orientation): 
	
	N = image.shape[0] 
	image_hat = np.fft.fft2(image)
	#print(image_hat.shape)
	image_hat = np.fft.fftshift(image_hat)
	#print(image_hat.shape)
	wx_j, wy_j = np.meshgrid(np.arange(-(N-1)/2, (N-1)/2 +1, dtype= int), np.arange(-(N-1)/2, (N-1)/2 +1,dtype= int))
	#image_hat = ((-1)**np.abs(wx_j+wy_j)/(N**2)*image_hat)
	
	image_hat = np.tile(image_hat[...,np.newaxis], [1,1,N])
	#print(image_hat.shape)
	#sys.exit()
	wz_j = np.arange(-(N-1)/2, (N-1)/2 +1, dtype= int)
	l_j = N * np.sinc(N * np.pi * wz_j)
	print("next")
	#l_j_hat = np.fft.fft(l_j)
	l_j_hat = np.tile(l_j[np.newaxis,np.newaxis,...],[N,N,1])
	print(l_j_hat.shape)
	#sys.exit()
	print("tiling of l_j done")
	#sys.exit()


	b_j_hat = image_hat * l_j_hat
	print("b_j_hat created")

	#sys.exit()
	b_j_hat = np.fft.ifftshift(b_j_hat)
	b_j = np.fft.ifftn(b_j_hat)


	#b_j =((-1)**np.abs(wz_j)*(N**3)*b_j)
	#print(b_j)
	#l_j = np.sin(wz_j * N * np.pi)/(np.pi * wz_j)
	return np.float32(np.real(b_j))


zika_file = mrcfile.open('zika_153.mrc')
rho = zika_file.data
R= RM.get_rotation_matrix()
image = project_fst(rho, R)
#print(image.shape)
print("Now getting image")
back_img = reconstruct(image, R)
print("Now saving")
with mrcfile.new('tmp.mrc', overwrite=True) as mrc:
	mrc.set_data(back_img)

#print(back_img.shape)
#with file('back_img.mrc', 'w') as outfile:
#    for slice_2d in image:
#        np.savetxt(outfile, slice_2d)
#np.savetxt('back_img.mrc', back_img, delimiter=',')
#back_obj = open("back_img.mrc", "w")
#back_obj.write(back_img)
#back_obj.close()


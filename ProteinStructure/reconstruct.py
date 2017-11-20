import numpy as np
from project_fst import project_fst
import mrcfile
import Rotation_Matrix as RM


def reconstruct(image, orientation): 
	
	N = image.shape[0] 
	image_hat = np.fft.fftn(image)
	image_hat = np.fft.fftshift(image_hat)

	wx_j, wy_j, wz_j = np.meshgrid(np.arange(-(N-1)/2, (N-1)/2 +1, dtype= int), np.arange(-(N-1)/2, (N-1)/2 +1,dtype= int),np.arange(-(N-1)/2, (N-1)/2 +1,dtype= int))
	image_hat = ((-1)**np.abs(wx_j+wy_j)/(N**2)*image_hat)


	#print(wz_j)
	l_j = N * np.sinc(N * np.pi * wz_j)
	l_j_hat = np.fft.fftn(l_j)
	#print(l_j_hat)
	b_j_hat = np.multiply(image_hat, l_j_hat)
	b_j_hat = np.fft.ifftshift(b_j_hat)
	b_j_hat = np.fft.ifftn(b_j_hat)


	b_j_hat =((-1)**np.abs(wx_j+wy_j)*(N**3)*b_j_hat)

	#print(np.sin(np.pi * wz_j * N))
	#l_j = np.sin(wz_j * N * np.pi)/(np.pi * wz_j)
	return b_j_hat
	#print("Done") 
	#l_j_hat = np.array() 

zika_file = mrcfile.open('zika_153.mrc')
rho = zika_file.data
R= RM.get_rotation_matrix()
image = project_fst(rho, R)

back_img = reconstruct(image, R)
#with file('back_img.mrc', 'w') as outfile:
#    for slice_2d in image:
#        np.savetxt(outfile, slice_2d)
#np.savetxt('back_img.mrc', back_img, delimiter=',')
#back_obj = open("back_img.mrc", "w")
#back_obj.write(back_img)
#back_obj.close()


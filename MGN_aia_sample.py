import matplotlib.pyplot as plt 
from scipy.io import readsav
import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter

dpath_SDO = '/Users/souvikb/Downloads/SDO_data_MGN/'
aia_171 = fits.getdata(dpath_SDO+'aia_l2_20210920_112027_3620108077_193.fits',ext=0) #This is the L2 SDO data obtained from co-aligned IRIS data
# The dimensions are in order (time, nx, ny)

#Writing the Multi-sale Gaussian Normalization (MGN) script for one AIA 171 image at a given time. 
dim_image = aia_171.shape
a0 = aia_171.min()
a1 = aia_171.max()
h=0.7
kernel_width = [1.25, 2.5, 5, 10, 20, 40]
Final_output = np.zeros((dim_image[0],dim_image[1],dim_image[2]))
for time in range(dim_image[0]):

	test_aia171 = aia_171[time,:,:] # image at 15th time step
	C_g = ((test_aia171 - a0)/(a1-a0))**(1/3.2) # Equation (4) of Morgan et al. (2014)
	C_prime = np.zeros((dim_image[1],dim_image[2],len(kernel_width)))

	for index in range(len(kernel_width)):
		loc_mean = gaussian_filter(test_aia171,sigma = kernel_width[index])
		diff_w = (test_aia171 - loc_mean)**2
		sigma_w = np.sqrt(gaussian_filter(diff_w,sigma = kernel_width[index]))
		C = (test_aia171 - loc_mean)/sigma_w
		#print(np.isnan(C.flatten()))
		C_prime[:,:,index] = np.arctan(0.7*C) # arctan(kC); where k =0.7 as per Morgan et al. (2014)

	Final_output[time,:,:] = h*C_g + ((1-h)/len(kernel_width))*np.nanmean(C_prime,axis=2)

	fig, ax = plt.subplots(1,2,figsize=(11,8))
	fig.subplots_adjust(wspace=0,left=0.1,right=0.9,top=0.9,bottom=0.1)
	ax=ax.ravel()

	ax[0].imshow(C_g, origin='lower')
	ax[0].set_title('Contrast Enhanced')
	ax[1].imshow(Final_output[time,:,:],origin='lower')
	ax[1].set_title('MGN Enhanced image')
	ax[1].xaxis.set_visible(False)
	ax[1].yaxis.set_visible(False)
	plt.show()
	plt.savefig(dpath_SDO+'Image_comparison-'+str(time)+'.png',dpi=400,bbox_inches = 'tight',pad_innches=0.1)

# plt.figure()
# plt.imshow(Final_output_193,origin='lower')
# plt.title('MGN enhanced')
# plt.show()

## rgb = np.dstack((Final_output_304,Final_output_171,Final_output_193)) # to make composite stacked images
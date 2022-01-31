import matplotlib.pyplot as plt 
from scipy.io import readsav
import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter
from matplotlib import colors
import cmasher as cmr
#import sunpy.cm as cm #The functionality of the sunpy.cm module is now in sunpy.visualization.colormaps as of SunPy 1.1
from tqdm import tqdm
from helita.io import lp
import COCOpy as cp
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.measure import label, regionprops
import matplotlib as mpl
from scipy import ndimage
import COCOpy as cp
import time
import multiprocessing as mp
import sunpy.cm as cm

dpath_SST = '/mn/stornext/d18/lapalma/reduc/2021/2021-08-04/CHROMIS/cubes_nb/'
dpath_SST_processed = '/mn/stornext/d18/lapalma/reduc/2021/2021-08-04/CHROMIS/cubes_nb/souvikb_CBP/'
dpath_EDVARDA_SDO = '/mn/stornext/d18/lapalma/reduc/2021/2021-08-04/CHROMIS/edvarda/sdo2sst/' 
dpath_SDO_BOSE ='/mn/stornext/d18/lapalma/reduc/2021/2021-08-04/CHROMIS/cubes_nb/souvikb_CBP/SDO_only_aia/'
sdo_target_cubes = '/mn/stornext/d18/lapalma/reduc/2021/2021-08-04/CHROMIS/edvarda/sdo/target/cubes/'

#----- time and wavelength stamps --------
time_steps = readsav(dpath_SST+'nb_4846_2021-08-04T09:56:50_scans=0-95_corrected_cmapcorr_im_times+wvl.idlsave')
times = time_steps['times']
wav = time_steps['LAMBDA']

hf1 = h5py.File(dpath_SST+'H_beta_widths_shifts_04.08.2021_CHROMIS.hdf5', 'r')
width = hf1['Width']
shift = hf1['Shift']
#cube_dopp=np.swapaxes(cube_dopp,0,1)

time_steps = readsav(dpath_SST+'nb_4846_2021-08-04T09:56:50_scans=0-95_corrected_cmapcorr_im_times+wvl.idlsave')
times = time_steps['times']
wav = time_steps['LAMBDA']
point1_aia = [400*0.037,600*0.037]
point2_aia = [1100*0.037,600*0.037]
point3_aia = [1100*0.037,1120*0.037]
point4_aia = [400*0.037,1120*0.037]

x_values_aia = [point1_aia[0], point2_aia[0]]
y_values_aia = [point1_aia[1], point2_aia[1]]

x_values1_aia =[point3_aia[0],point4_aia[0]]
y_values1_aia =[point3_aia[1],point4_aia[1]]

x_values2_aia =[point1_aia[0],point1_aia[0]]
y_values2_aia =[point1_aia[1],point3_aia[1]]

x_values3_aia =[point2_aia[0],point2_aia[0]]
y_values3_aia =[point3_aia[1],point2_aia[1]]

for time_index in range(96):
	fig, axs =plt.subplots(1,2,figsize=(12,7.5),facecolor='w', edgecolor='k')
	fig.subplots_adjust(hspace = 0.0,wspace=0.08,left=0.1,right=0.9,top=0.9,bottom=0.1)
	axs=axs.ravel()

	imc1=axs[0].imshow(width[:,:,time_index],origin='lower',vmax=0.82,vmin=0.45,cmap='sdoaia304',extent=[0,1796*0.037,0,1138*0.037])
	dividerc1 = make_axes_locatable(axs[0])
	caxc1 = dividerc1.append_axes("top", size="5%", pad=0.05)
	cbar2=plt.colorbar(imc1,cax=caxc1,orientation='horizontal',)
	#cbar2.set_label(r'H$\beta$ LC width [$\AA$]',size=12,color='black',weight='bold')
	axs[0].text(.05, .9, r'H$\beta$ LC width [$\AA$]', transform=axs[0].transAxes, fontsize=13,color='black',weight='bold')
	cbar2.ax.xaxis.set_ticks_position("top")
	axs[0].set_xlabel('X$_{\mathrm{1}}$ [arcsec]')
	axs[0].set_ylabel('Y$_{\mathrm{1}}$ [arcsec]')
	axs[0].set_yticks([0,10,20,30,40])
	b = times[time_index][0:8]
	axs[0].text(.75, .06, b.decode('utf-8')+' UTC',ha="center",color='black',size=12,transform=axs[0].transAxes,weight='bold')
	axs[0].plot(x_values_aia, y_values_aia,color='white',linestyle='dashed')
	axs[0].plot(x_values1_aia, y_values1_aia,color='white',linestyle='dashed')
	axs[0].plot(x_values2_aia, y_values2_aia,color='white',linestyle='dashed')
	axs[0].plot(x_values3_aia, y_values3_aia,color='white',linestyle='dashed')



	imc2=axs[1].imshow((shift[:,:,time_index])*62,origin='lower',vmax=6,vmin=-6,cmap='RdBu_r',extent=[0,1796*0.037,0,1138*0.037])
	dividerc1 = make_axes_locatable(axs[1])
	caxc1 = dividerc1.append_axes("top", size="5%", pad=0.05)
	cbar2=plt.colorbar(imc2,cax=caxc1,orientation='horizontal',)
	#cbar2.set_label(r'V$_{\mathrm{Dopp}}$ [km s^$^{-1}$]',size=12,color='black',weight='bold')
	axs[1].text(.05, .9, r'V$_{\mathrm{LOS}}$ [km s$^{-1}$]', transform=axs[1].transAxes, fontsize=13,color='black',weight='bold')
	cbar2.ax.xaxis.set_ticks_position("top")
	axs[1].set_xlabel('X$_{\mathrm{1}}$ [arcsec]')
	axs[1].set_yticks([0,10,20,30,40])
	axs[1].plot(x_values_aia, y_values_aia,color='black',linestyle='dashed')
	axs[1].plot(x_values1_aia, y_values1_aia,color='black',linestyle='dashed')
	axs[1].plot(x_values2_aia, y_values2_aia,color='black',linestyle='dashed')
	axs[1].plot(x_values3_aia, y_values3_aia,color='black',linestyle='dashed')

	plt.savefig('/mn/stornext/d9/souvikb/CBP_spicules_SDO/widths_shifts/lcwidths_shifts-'+str(time_index)+'.png',dpi=300,bbox_inches = 'tight',pad_innches=0.1)
	#plt.show()
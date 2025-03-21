#*****************************************************************************
#                       COLD MASK IMAGE PROGRAM   [Rix+2004]
#-----------------------------------------------------------------------------
#                             Cristina Martinez 08/20
#*****************************************************************************

# Versions
# %pip install numpy== 1.18.5
# %pip install matplotlib==3.2.2
# %pip install astropy==4.0
# %pip install photutils==0.7.2

#!pip show numpy
#!pip show matplotlib
#!pip show astropy
#!pip show photutils


# Basics
import numpy as np

# Matplotlib
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
plt.close('all')
plt.ion()

# Astropy
from photutils import detect_threshold
from astropy.convolution import Gaussian2DKernel
from photutils import detect_sources
from photutils import deblend_sources
from astropy.convolution import convolve

# CML
import sys
#sys.path.append('/Users/felipe/Dropbox/CML/')
sys.path.append('/Users/z3530031/Dropbox/CML/')

from CML import DeleteMaskRegion
from CML import Patchs



# =============================================================================
# Esta lineas es para saber el tiempo de ejecución
# =============================================================================

import time
start = time.time()


# =============================================================================
# Load images and masks -- EDIT AS REQUIRED
# =============================================================================

# Group info
groupID = '400138'
band = 'gri'  


# Paths & data
masks_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/masks/intermediate/'
cold = fits.getdata(masks_path+'ColdMasked.fits') 

gri = fits.getdata('/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/data/images/'+groupID+'-gri.fits')
header = fits.getheader('/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/data/images/'+groupID+'-gri.fits')



# =============================================================================
# Unsharped image
# =============================================================================
# Convolve the image with a gaussian kernel to create an unsharp image
sigma_pix = 5
kernel= Gaussian2DKernel(sigma_pix) #X & Y size by default 8*sigma_pix
conv_gri = convolve(gri, kernel)

# Create cold mask
cold_mask = np.zeros_like(cold)
cold_mask[np.isnan(cold)] = np.nan
cold_mask = np.isnan(cold_mask)

# Subtract data - conv_data
unsharp = gri - conv_gri

# Masking with the Cold mask
data = np.zeros_like(unsharp)
data[cold_mask] = np.nan
data[~cold_mask] = unsharp[~cold_mask]


# =============================================================================
# Sorce detection & deblending en secciones mas pequenas sobre la imagen
# =============================================================================
patchs = Patchs(data,1800)

# Construimos una matriz de 0 del tamaño de la imagen como estructura para la mascara hot
mask = np.zeros(np.shape(data))

for i in patchs:
	# Establish a Hot nsigma threshold (we detect everything above)
	threshold = detect_threshold(data[i], nsigma=1.2, sigclip_sigma=3., mask_value = np.nan)
	# Segmentation - detecting sources **********************************************
	segm = detect_sources(data[i], threshold, npixels=7, mask=cold_mask[i])
	# Deblending sources **********************************************
	#mask[i] = deblend_sources(data[i], segm, npixels=7, filter_kernel=kernel).data
	mask[i] = deblend_sources(data[i], segm, npixels=7).data


# =============================================================================
# Fin contador de tiempo
# =============================================================================
print('It took', time.time()-start, 'seconds.')


# =============================================================================
# Remove masks over the target if required
# =============================================================================
print('Remove masks over the BCG')

patchs = Patchs(data,1800) # Use a larger region than before to see the 'borders'


for i in patchs:
    #Delete 'by hand' clicking over the required position and we replace in the original image
    mask[i] = DeleteMaskRegion(mask[i],data[i], vmax=1.5)




# =============================================================================
# Convolve the mask with a Gaussian kernel to make it larger 
# =============================================================================

sigma_pix = 1
kernel= Gaussian2DKernel(sigma_pix, x_size=3, y_size=3) 

conv_mask = convolve(mask, kernel)


# =============================================================================
# Mask the image with the final mask
# =============================================================================

data[conv_mask != 0]= np.nan
cold[conv_mask != 0] = np.nan



# =============================================================================
# Save masked image
# =============================================================================

masked_hdu = fits.PrimaryHDU(cold)  # we create a PrimaryHDU object to encapsulate the data:
wcs = WCS(header)
masked_hdu.header.update(wcs.to_header())
masked_hdu.header.append(('HISTORY', 'Hot + Cold Mask', 'C. Martinez-Lombilla'), end=True)
masked_hdu.writeto(masks_path+'HotColdMasked.fits', overwrite=True) # write to a new file



# =============================================================================
# Save unshaped masked image
# =============================================================================


masked_hdu = fits.PrimaryHDU(data)  # we create a PrimaryHDU object to encapsulate the data:
wcs = WCS(header)
masked_hdu.header.update(wcs.to_header())
masked_hdu.header.append(('HISTORY', 'Hot + Cold Mask Ushaped', 'C. Martinez-Lombilla'), end=True)
masked_hdu.writeto(masks_path+'UnshapedMasked.fits', overwrite=True) # write to a new file

















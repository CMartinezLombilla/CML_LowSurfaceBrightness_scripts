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

import sys

# CML
sys.path.append('/Users/z3530031/Dropbox/CML/')

from CML import DeleteMaskRegion
from CML import Patchs



# =============================================================================
# Esta lineas es para saber el tiempo de ejecución
# =============================================================================

import time
start = time.time()


# =============================================================================
# Cargamos las imagenes y las mascaras 
# =============================================================================

# Group info
groupID = '400138'
band = 'gri'  

# Paths & data
masks_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/masks/intermediate/'

image = fits.getdata('/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/data/images/'+groupID+'-gri.fits')
header = fits.getheader('/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/data/images/'+groupID+'-gri.fits')

cold_input = fits.getdata('masks/cold_input.fits')

#================================================================================
# We mask the data with the input mask from the bright star ds9 regions
#================================================================================

data = np.zeros_like(image)
data[cold_input!=0] = np.nan
data[cold_input==0] = image[cold_input==0]


#================================================================================
# COLD Source detection & deblending -- bright and extended sources
#================================================================================

# Gaussian kernel for convolve the data for a better source detection
sigma_pix = 5
kernel = Gaussian2DKernel(sigma_pix) #X & Y size by default 8*sigma_pix
kernel.normalize()

# =============================================================================
# Define smaller section over the image it's is too big -- image size if not needed
# =============================================================================
patchs = Patchs(data,1800)

mask = np.zeros(np.shape(data))

for i in patchs:
    # Establish a Hot nsigma threshold (we detect everything above)
    threshold = detect_threshold(data[i], nsigma=1.1, sigclip_sigma=3., mask_value = np.nan)

    # Segmentation - detecting sources
    segm = detect_sources(data[i], threshold, npixels=40,filter_kernel=kernel)
    # Deblending sources 
    mask[i] = deblend_sources(data[i], segm, npixels=40, filter_kernel=kernel).data





# =============================================================================
# Fin contador de tiempo
# =============================================================================
print('It took', time.time()-start, 'seconds.')


# =============================================================================
# Remove mask we don't need
# =============================================================================
print('Remove masks over the BCG')

patchs = Patchs(data,1800) # Use a larger region than before to see the 'borders'

for i in patchs:
    # Editamos a mano la seccion y la reemplazamos en la original
    mask[i] = DeleteMaskRegion(mask[i],data[i], vmax=0.5)



# =============================================================================
# Convolve the mask with a Gaussian kernel to make it larger 
# =============================================================================

sigma_pix = 2   
kernel= Gaussian2DKernel(sigma_pix)  #X & Y size by default 8*sigma_pix

conv_mask = convolve(mask, kernel)

# =============================================================================
# Mask the image with the final mask
# =============================================================================

data[conv_mask != 0]= np.nan


# =============================================================================
# Save masked image
# =============================================================================

# Create a fits file to save the  mask 
masked_hdu = fits.PrimaryHDU(data)  # we create a PrimaryHDU object to encapsulate the data:
wcs = WCS(header)
masked_hdu.header.update(wcs.to_header())
masked_hdu.header.append(('HISTORY', 'Cold Mask', 'C. Martinez-Lombilla'), end=True)
masked_hdu.writeto(masks_path+'ColdMasked.fits', overwrite=True) # write to a new file



#*****************************************************************************
#                       RETOUCH MASK IMAGE PROGRAM   [Rix+2004]
#-----------------------------------------------------------------------------
#                             Cristina Martinez 04/21
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

# Las siguiente dos lineas son para que no pete el MATPLOTLIB fuera de ANACONDA
import matplotlib
matplotlib.use('TkAgg')

# Matplotlib
import matplotlib.pyplot as plt
plt.close('all')
plt.ion()

# Astropy
from astropy.io import fits
from astropy.wcs import WCS


# CML
import sys
sys.path.append('/Users/z3530031/Dropbox/CML/')

from CML import CircularHandMask
from CML import Patchs



# =============================================================================
# Load images and masks -- EDIT AS REQUIRED
# =============================================================================

# Group info
groupID = '400138'
band = 'gri'  


# Paths & data
masks_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/masks/intermediate/'
unsh = fits.getdata(masks_path+'UnshapedMasked.fits')
unsh_head = fits.getheader(masks_path+'UnshapedMasked.fits')

file = masks_path+'UnshapedMasked'


# =============================================================================
# retouch unshaped image by adding masks
# =============================================================================
print('Unshaped: Add manual masks of a given size over NEGATIVE regions around sources')

patchs = Patchs(unsh,600) # Change size as necessary -- or equal to image size if the image is small

for i in patchs:
    unsh[i] = CircularHandMask(unsh[i], mask_value= np.nan, vmin = 0., vmax=0.3)



# =============================================================================
# Save masked image
# =============================================================================

# Create a fits file to save the  mask 
masked_hdu = fits.PrimaryHDU(unsh)  # we create a PrimaryHDU object to encapsulate the data:
wcs = WCS(unsh_head)
masked_hdu.header.update(wcs.to_header())
masked_hdu.writeto(file+'-retouched.fits', overwrite=True) # write to a new file

















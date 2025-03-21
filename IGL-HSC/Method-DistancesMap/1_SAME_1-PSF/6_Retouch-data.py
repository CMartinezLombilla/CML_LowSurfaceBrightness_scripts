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
sys.path.append('/Users/felipe/Dropbox/CML/')
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
masks_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/masks/'
data = fits.getdata(masks_path+'intermediate/HotColdMasked.fits')
data_head = fits.getheader(masks_path+'intermediate/HotColdMasked.fits')

file = masks_path+'intermediate/HotColdMasked'


# =============================================================================
# First iteration -- comment if doing 'following iterations'

unsh = fits.getdata(masks_path+'intermediate/UnshapedMasked-retouched.fits')
# Mask loaded data with Hot + Cold + retouched (unsh) masks
data[np.isnan(unsh)] = np.nan
# =============================================================================


'''
# =============================================================================
# Following iterations - uncomment if want to retouch againg the retouched image
 
retou = fits.getdata(masks_path+'intermediate/HotColdMasked-retouched.fits')
# Mask loaded data with Hot + Cold + retouched (unsh) masks
data[np.isnan(retou)] = np.nan
# =============================================================================
'''


# =============================================================================
# Retouch image by adding masks
# =============================================================================
print('gri image: Add manual masks of a given size over BRIGHT regions around sources')

patchs = Patchs(data,2800)  # Change size as necessary -- or equal to image size if the image is small

for i in patchs:
    data[i] = CircularHandMask(data[i], mask_value= np.nan, vmin = 0., vmax=0.5) # Change vmin & max as necessary 



# =============================================================================
# Save masked image
# =============================================================================

# Create a fits file to save the  mask 
masked_hdu = fits.PrimaryHDU(data)  # we create a PrimaryHDU object to encapsulate the data:
wcs = WCS(data_head)
masked_hdu.header.update(wcs.to_header())
masked_hdu.writeto(file+'-retouched.fits', overwrite=True) # write to a new file



# =============================================================================
# Create and save a FINAL mask -- comment if want to keep retouching
# =============================================================================

final_mask = np.zeros_like(data)
final_mask[np.isnan(data)] = np.nan

mask_hdu = fits.PrimaryHDU(final_mask) 
wcs = WCS(data_head)
mask_hdu.header.update(wcs.to_header())
mask_hdu.header.append(('HISTORY', ' Final BCG Mask', 'C. Martinez-Lombilla'), end=True)
mask_hdu.writeto('masks/PSFMask.fits', overwrite=True) # write to a new file


masked_hdu = fits.PrimaryHDU(data)  
wcs = WCS(data_head)
masked_hdu.header.update(wcs.to_header())
masked_hdu.writeto('masks/PSFMasked.fits', overwrite=True) # write to a new file

















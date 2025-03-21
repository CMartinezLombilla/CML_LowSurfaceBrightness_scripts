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

from CML import HandMask
from CML import Patchs



# =============================================================================
# Cargamos las imagenes y las mascaras 
# =============================================================================

# Group info
groupID = '400138'
band = 'gri'  


# Paths & data
masks_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/masks/intermediate/'
unsh = fits.getdata(masks_path+'UnshapedMasked-retouched.fits')
unsh_head = fits.getheader(masks_path+'UnshapedMasked.fits')

file = masks_path+'UnshapedMasked'


# =============================================================================
# Retocamos la imagen unshaped agregando mascaras a mano
# =============================================================================
print('Unshaped: Add manual masks of a given size over NEGATIVE regions around sources')

patchs = Patchs(unsh,1500)

for i in patchs:
    unsh[i] = HandMask(unsh[i], mask_value= np.nan, vmin = 0., vmax=0.3)



# =============================================================================
# Guardamos la nueva imagen enmascarada 
# =============================================================================

# Create a fits file to save the  mask 
masked_hdu = fits.PrimaryHDU(unsh)  # we create a PrimaryHDU object to encapsulate the data:
wcs = WCS(unsh_head)
masked_hdu.header.update(wcs.to_header())
masked_hdu.writeto(file+'-retouched.fits', overwrite=True) # write to a new file
















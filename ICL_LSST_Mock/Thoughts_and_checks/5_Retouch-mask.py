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
import matplotlib.pyplot as plt
plt.close('all')
plt.ion()

# Astropy
from astropy.io import fits
import os
import glob


# CML
from CML import CircularHandMask
from CML import Patchs
from CML import save_fits_image



# =============================================================================
# Initial params
# =============================================================================
sim = 'Horizon_AGN'
cluster = '048'
orientation = 'xy'
iteration = 1

# =============================================================================
# Paths 
# =============================================================================
mask_dir = '/Users/z3530031/unsw/Work/Projects/LSST/Chall3-ICL/results/'+sim+'/'+cluster+'_'+orientation+'/Masks/'

new_file_name = 'retouched'

'''
# =============================================================================
# First iteration -- comment if doing 'following iterations'

im_files = glob.glob(mask_dir+'combined*.fits')
im, hdr = fits.getdata(im_files[0], header = True) 
# =============================================================================
'''

# =============================================================================
# Following iterations - uncomment if want to retouch againg the retouched image

iteration = 2
 
im_files = glob.glob(mask_dir+'retouch'+str(iteration-1)+'*.fits')
im, hdr = fits.getdata(im_files[0], header = True) 
# =============================================================================



# =============================================================================
# Retouch image by adding masks
# =============================================================================
print('Image: Add manual masks of a given size over BRIGHT regions around sources')

patchs = Patchs(im,1500)  # Change size as necessary -- or equal to image size if the image is small

for patch in patchs:
    im[patch] = CircularHandMask(im[patch], mask_value=np.nan, vmin=0., vmax=20) # Change vmin & max as necessary 


# =============================================================================
# Save new masked.fits file
# =============================================================================
save_fits_image(im, hdr, mask_dir, 'retouch'+str(iteration)+'-'+os.path.basename(im_files[0]))
    

'''
# =============================================================================
# Create and save a FINAL mask -- comment if want to keep retouching
# =============================================================================
save_fits_image(im, hdr, mask_dir, 'final-'+sim+'-'+cluster+'-'+orientation+'.fits')
'''














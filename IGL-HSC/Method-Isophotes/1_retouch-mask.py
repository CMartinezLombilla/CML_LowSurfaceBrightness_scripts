#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 15:43:36 2022

@author: C. Martinez-Lombilla
"""

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


# CML
import glob
import sys
sys.path.append('/Users/felipe/Dropbox/CML/')
sys.path.append('/Users/z3530031/Dropbox/CML/')

from CML import CircularHandMask_new
from CML import Patchs
from CML import save_fits_image


# =============================================================================
# Load images and masks -- EDIT AS REQUIRED
# =============================================================================

# Group info
groupID = '400138'
# Bands
bands = ['g', 'r', 'i']  


# Paths & data
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'
masks_path = results_path+'masks/'

data_iso = fits.getdata(results_path+'Isophotes/isophotes_i.fits')


# =============================================================================
# Retouch image by adding masks
# =============================================================================
print('Add manual masks of a given size over BRIGHT regions around sources')

patchs = Patchs(data_iso,900)  # Change size as necessary -- or equal to image size if the image is small
new_mask = np.zeros_like(data_iso) 

for i in patchs:
    new_mask[i] = CircularHandMask_new(data_iso[i], mask_value= np.nan, vmin = 5., vmax=32) # Change vmin & max as necessary 

 

# =============================================================================
# Apply mask over all files
# =============================================================================

image_paths = glob.glob(masks_path+'*.fits')

for im_path in image_paths:

    # Load images
    data, hdr = fits.getdata(im_path, header = True)
    
    data[np.isnan(new_mask)] = np.nan

    save_fits_image(data, hdr, masks_path+'test/', im_path[82:])
    
    












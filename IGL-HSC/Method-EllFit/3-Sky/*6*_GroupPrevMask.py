#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 15:30:11 2020

@author: C. Martinez-Lombilla
"""

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
import matplotlib.pyplot as plt
plt.close('all')
plt.ion()

# Astropy
from astropy.io import fits
from astropy.wcs import WCS

# CML
import sys
#sys.path.append('/Users/felipe/Dropbox/CML/')
sys.path.append('/Users/z3530031/Dropbox/CML/')

from CML import HandMask
from CML import Patchs

# =============================================================================
# Cargamos las imagenes y las mascaras 
# =============================================================================
# Group info
groupID = '400282'
band = 'gri'  

main_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/HSC-Mask/'+groupID
file = main_path+'/0_data/gri-cut'


# Cargamos los datos
data = fits.getdata(file+'.fits') 
groupMask = np.zeros_like(data)



# =============================================================================
# Ponemos mascaras a mano sobre las galaxias del cumulo
# =============================================================================
patchs = Patchs(data,2900) 

for i in patchs:
	data[i] = HandMask(data[i], mask_value= np.nan, vmax=0.5)

groupMask[np.where(np.isnan(data))] = np.nan


# =============================================================================
# Guardamos la imagen enmascarada 
# =============================================================================

# Create a fits file to save the  mask 
groupMask_hdu = fits.PrimaryHDU(groupMask)  # we create a PrimaryHDU object to encapsulate the data:
wcs = WCS(fits.getheader(file+'.fits'))
groupMask_hdu.header.update(wcs.to_header())
groupMask_hdu.header.append(('HISTORY', 'Group Mask', 'C. Martinez-Lombilla'), end=True)
groupMask_hdu.writeto(main_path+'/0_masks/intermediate/GroupOnlyMask.fits', overwrite=True) # write to a new file



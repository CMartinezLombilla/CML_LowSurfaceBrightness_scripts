#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 16:45:25 2020

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
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from photutils import detect_sources
from photutils import deblend_sources
from photutils import detect_threshold

# CML
import sys
#sys.path.append('/Users/felipe/Dropbox/CML/')
sys.path.append('/Users/z3530031/Dropbox/CML/')

from CML import HandMask
from CML import Patchs
from CML import DeleteMaskRegion


# =============================================================================
# Cargamos las imagenes y las mascaras 
# =============================================================================
# =============================================================================
# Cargamos las imagenes y las mascaras 
# =============================================================================
# Group info
groupID = '400282'
band = 'gri'  

main_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/HSC-Mask/'+groupID
file = main_path+'/0_data/gri-cut'

# Cargamos los datos
gri =  fits.getdata(file+'.fits') 
BCGmask = fits.getdata(main_path+'/0_masks/BCGMask.fits') 
groupOnlyMask = fits.getdata(main_path+'/0_masks/intermediate/GroupOnlyMask.fits')

# Creamos mascara hot&cold - parches sobre las galaxias del grupo
group_mask = np.empty_like(gri)

group_mask[np.isnan(BCGmask)] = np.nan
group_mask[np.isnan(groupOnlyMask)]= 0
group_mask[~np.isnan(group_mask)]= 0
group_mask = np.isnan(group_mask)


# Convolucionamos imagen gri con kernel gaussiano [to creat an unsharp image]
sigma_pix = 5
kernel= Gaussian2DKernel(sigma_pix) #X & Y size by default 8*sigma_pix
conv_gri = convolve(gri, kernel)

# Creamos unsharped image: restamos gri - conv_gri
unsharp = gri - conv_gri

# enmascaramos unsharped con la mascara del grupo
data = np.zeros_like(unsharp)
data[group_mask] = np.nan
data[~group_mask] = unsharp[~group_mask]


# =============================================================================
# Sorce detection & deblending para detectarfuentes alrededor de las galaxias 
# =============================================================================
patchs = Patchs(data,600)

# Construimos una matriz de 0 del tamaño de la imagen como estructura para la mascara hot
mask = np.zeros(np.shape(data))

for i in patchs:
	# Establish a Hot nsigma threshold (we detect everything above)
	threshold = detect_threshold(data[i], nsigma=1.2, sigclip_sigma=3., mask_value = np.nan)
	# Segmentation - detecting sources **********************************************
	segm = detect_sources(data[i], threshold, npixels=7, mask=group_mask[i])
	# Deblending sources **********************************************
	#mask[i] = deblend_sources(data[i], segm, npixels=7, filter_kernel=kernel).data
	mask[i] = deblend_sources(data[i], segm, npixels=7).data



# =============================================================================
# Quitamos máscaras de sobra
# =============================================================================
print('Remove masks over the galaxies')

patchs = Patchs(data,870) # Cambiamos el tamaño de los parches


for i in patchs:
    # Editamos a mano la seccion y la reemplazamos en la original
    mask[i] = DeleteMaskRegion(mask[i],data[i], vmax=1.5)



# =============================================================================
# Retocamos la imagen agregando mascaras a mano
# =============================================================================
data[mask != 0]= np.nan

data_masked = gri
data_masked[np.isnan(data)] = np.nan

patchs = Patchs(data,870)

for i in patchs:
    data_masked[i] = HandMask(data_masked[i], mask_value= np.nan, vmax=10)



# =============================================================================
# Guardamos mascara e imagen enmascarada 
# =============================================================================

masked_hdu = fits.PrimaryHDU(data_masked)  # we create a PrimaryHDU object to encapsulate the data:
wcs = WCS(fits.getheader(file+'.fits'))
masked_hdu.header.update(wcs.to_header())
masked_hdu.writeto(main_path+'/0_masks/GroupMasked.fits', overwrite=True) # write to a new file


group_mask2 = np.zeros_like(data_masked) 
group_mask2[np.isnan(data_masked)]= np.nan

mask_hdu = fits.PrimaryHDU(group_mask2)  # we create a PrimaryHDU object to encapsulate the data:
wcs = WCS(fits.getheader(file+'.fits'))
mask_hdu.header.update(wcs.to_header())
mask_hdu.writeto(main_path+'/0_masks/GroupMask.fits', overwrite=True) # write to a new file

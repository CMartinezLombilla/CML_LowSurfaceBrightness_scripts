#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 11:59:44 2021

@author: C. Martinez-Lombilla

=============================================================================
    Programa que genera combina todas las mascaras bool de cada region 
    (asociada a una bright star). Input para 'cold' mask.
=============================================================================
"""


import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import glob

# =============================================================================
# Leemos las mascaras de todas las regiones y cargamos sus datos
# =============================================================================
results_dir = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'

masks = glob.glob(results_dir+'PSF/masks/mask_s*.fits')
data = []
header = []

for star_file in masks:
    data.append(fits.getdata(star_file))
    header.append(fits.getheader(star_file))

# =============================================================================
# Creamos la máscara combinada [suma!]
# =============================================================================

new_mask = np.zeros_like(data[0])
for star_data in data:
    new_mask += star_data


mask_hdu = fits.PrimaryHDU(new_mask)  # we create a PrimaryHDU object to encapsulate the data
wcs = WCS(header[0])
mask_hdu.header.update(wcs.to_header())
mask_hdu.writeto(results_dir+'PSF/masks/cold_input.fits', overwrite=True) # write to a new file


plt.imshow(new_mask,origin='lower')
plt.show()
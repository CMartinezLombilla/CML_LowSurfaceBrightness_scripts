#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 17:36:54 2021

@author: C. Martinez-Lombilla
"""

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from regions import read_ds9


# Group info
groupID = '400138'
bands = ['g', 'r', 'i']  

# Path
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'




for band in bands:

    # =============================================================================
    # Leemos las regiones a desenmascarar
    # =============================================================================
    reg = read_ds9(results_path+'masks/ds9Regions/StarsCorr.reg')
    OrigMask = fits.getdata(results_path+'masks/BCGMask-'+band+'.fits')
    header = fits.getheader(results_path+'masks/BCGMask-'+band+'.fits')
    
    # =============================================================================
    # Hacemos un bucle que pasa por las regiones para crear la máscara
    # =============================================================================
    
    mask = np.zeros_like(OrigMask)
    shape = mask.shape
    wcs = WCS(header)
    
    for i in reg:    
       mask = mask +i.to_pixel(wcs).to_mask('center').to_image(shape)
       
    
    # =============================================================================
    # Enmascaramos
    # =============================================================================
    OrigMask[mask!=0]=np.nan


    # =============================================================================
    # Guardamos la nueva imagen enmascarada 
    # =============================================================================
    
    # Create a fits file to save the  mask 
    masked_hdu = fits.PrimaryHDU(OrigMask)  # we create a PrimaryHDU object to encapsulate the data:
    wcs = WCS(header)
    masked_hdu.header.update(wcs.to_header())
    masked_hdu.writeto(results_path+'masks/BCGMask-'+band+'2.fits', overwrite=True) # write to a new file
    

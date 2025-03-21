#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 10 17:15:00 2021

@author: C. Martinez-Lombilla
"""


from astropy.io import fits
from astropy.wcs import WCS
import numpy as np


# =============================================================================
# Cargamos las imagenes y las mascaras 
# =============================================================================

# Group info
groupID = '400138'
bands = ['g', 'r', 'i'] 


# Path
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'
masks_path = results_path+'masks/'

for band in bands:
    # Data
    data = fits.getdata(results_path+'2D-Imfit/'+band+'/2DresSersPSF-'+band+'.fits')
    header = fits.getheader(masks_path+'GroupMasked-'+band+'.fits')

    wcs = WCS(header)
    
    new_hdu = fits.PrimaryHDU(data)  # we create a PrimaryHDU object to encapsulate the data
    new_hdu.header.update(wcs.to_header())
    new_hdu.writeto(results_path+'Dist_IGLlight-'+band+'.fits', overwrite=True) # write to a new file


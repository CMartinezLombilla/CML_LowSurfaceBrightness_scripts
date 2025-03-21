#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 10 17:15:00 2021

@author: C. Martinez-Lombilla
"""


from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

from CML import save_fits_image


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

    igl = fits.getdata(results_path+'2D-Imfit/'+band+'/distMap/2DresSersPSF-'+band+'.fits')
    header = fits.getheader(masks_path+'GroupMasked-'+band+'.fits')
    
    save_fits_image(igl, header, results_path, 'Dist_IGLlight-'+band+'.fits')
    
    bcg = fits.getdata(results_path+'2D-Imfit/'+band+'/distMap/BCGfrac-2DresPSF-'+band+'.fits')
    save_fits_image(bcg, header, results_path, 'BCG_and_IGLlight-'+band+'.fits')

        

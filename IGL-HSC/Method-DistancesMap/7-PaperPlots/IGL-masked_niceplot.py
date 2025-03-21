#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 16:29:05 2022

@author: C. Martinez-Lombilla
"""



from astropy.io import fits
from regions import read_ds9
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
          
    # Data
    igl = fits.getdata(results_path+'Dist_IGLlight-'+band+'.fits')
    Cmask = fits.getdata(masks_path+'CoreMask-'+band+'.fits')    
    
    header = fits.getheader(results_path+'Dist_IGLlight-'+band+'.fits')
    
    # Mask data
    igl_masked = np.copy(igl)
    igl_masked[np.isnan(Cmask)]=np.nan   

    
    # =============================================================================
    # Save IGL data masked
    # =============================================================================    
    # IGL
    save_fits_image(igl_masked, header, results_path, 'prety_IGLMasked-'+band+'.fits')

    

    

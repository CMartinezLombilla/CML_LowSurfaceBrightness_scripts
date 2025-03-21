#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 15:33:09 2021

@author: C. Martinez-Lombilla
"""

import os
import numpy as np


from astropy.io import fits
from astropy.wcs import WCS


# =============================================================================
# Parametros iniciales del grupo/imagen/IGL/aperturas
# =============================================================================
groupID = '400138'
bands = ['g', 'r', 'i']


# =============================================================================
# Paths
# =============================================================================
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'



# =============================================================================
# Diccionarios para cada imagen
# =============================================================================


files = {'g':results_path+'IGLlight-g.fits',
         'r':results_path+'IGLlight-r.fits', 
         'i':results_path+'IGLlight-i.fits'}



for band in bands:    
    
    data = fits.getdata(files[band])
    HSC_mask = fits.getdata(results_path+'masks/intermediate/detectorHSC-'+band+'.fits')

    data[HSC_mask>0] = np.nan   
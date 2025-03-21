#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 15:57:43 2022

@author: C. Martinez-Lombilla
"""

import numpy as np
import os
from astropy.io import fits
from astropy.nddata import block_reduce
import fnmatch

from CML import save_fits_image

# =============================================================================
# Initial params
# =============================================================================
sim = 'Horizon_AGN'
cluster = '048'
orientation = 'xy'
bands = ['r']


# =============================================================================
# Paths 
# =============================================================================
mask_dir = '/Users/z3530031/unsw/Work/Projects/LSST/Chall3-ICL/results/'+sim+'/'+cluster+'_'+orientation+'/Masks/'


# =============================================================================
# Load masks cold and unsharped
# =============================================================================
cold_file = fnmatch.filter(os.listdir(mask_dir), 'cold*.fits')[0]
unsharp_file = fnmatch.filter(os.listdir(mask_dir), 'unsharp*.fits')[0]


cold, hdr = fits.getdata(mask_dir+cold_file, header = True) 
unsharp = fits.getdata(mask_dir+unsharp_file) 



# =============================================================================
# Combine masks 
# =============================================================================
combined = np.copy(unsharp)
combined[np.isnan(cold)] = np.nan


# =============================================================================
# Save new .fits file
# =============================================================================
save_fits_image(combined, hdr, mask_dir, 'combined'+cold_file[4:])
    
    
    
    
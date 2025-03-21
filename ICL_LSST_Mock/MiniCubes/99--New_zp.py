#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 15:57:43 2022

@author: C. Martinez-Lombilla


!!!!!!!!!!!!!!!!!!!!!!!


CODE NOT FINISHED!!


!!!!!!!!!!!!!!!!!!!!!!!

"""

import numpy as np
import os
import glob

import matplotlib.pyplot as plt
plt.close('all')

from astropy.io import fits
from astropy.nddata import block_reduce
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

from CML import save_fits_image

# =============================================================================
# Initial params
# =============================================================================
sims = ['Horizon_AGN', 'Hydrangea', 'Magneticum', 'TNG_100'] 
orientation = 'xy'
bands = ['r']

npix_bin = 2 #number of pixels for the binning
factor = 1e12 # scale factor for multiplying the data values
factor_str = format(factor, '.0e')


for sim in sims:
    
    # =============================================================================
    # Paths 
    # =============================================================================
    data_dir = '/Users/z3530031/unsw/Work/Projects/LSST/Chall3-ICL/data/'+sim+'/'
    res_dir = '/Users/z3530031/unsw/Work/Projects/LSST/Chall3-ICL/results/'+sim+'/'
    
    #image_file = '00761_0000048_0.05_xy_r_4Mpc'
    
    raw_im_files = glob.glob(data_dir+'*xy_r_4Mpc.fits') # Load a list 
    
 
    for clust_files in raw_im_files:
        
        cluster = clust_files[-23:-20]
        
        
        for ind, band in enumerate(bands):
        
            # Load images
            new_im_dir = res_dir+cluster+'_'+orientation+'/Images/'
            new_filename = sim+'_'+cluster+'_'+orientation+'_'+band+'_2Mpc'+'_bin'+str(npix_bin)+'_x'+factor_str+'.fits'
            
            # New binned, xfactor, cut image
            new_im, new_hdr = fits.getdata(new_im_dir+new_filename, header = True)
            
            # Raw/original image
            raw_im, raw_hdr = fits.getdata(clust_files, header = True)
            
            # =============================================================================
            # OBTAIN ZEROPOINT
            # =============================================================================
            #raw_ct = counts in a given pixel or region in the raw image 
            #new_ct = counts in THE SAME pixel or region in the new image
            
            #zp = -2.5*np.log10(raw_ct/new_ct)
            
            # =============================================================================
            # SAVE ZEROPOINT VALUE IN THE HEADER   
            # =============================================================================


            #galaxy_fits=fits.open("rotated.fits")

            #new_image[0].header['HISTORY'] = 'Cut image'        
            
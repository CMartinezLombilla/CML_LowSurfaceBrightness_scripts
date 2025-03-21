#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 15:57:43 2022

@author: C. Martinez-Lombilla
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
sims = ['Hydrangea'] #'Horizon_AGN', 'Hydrangea', 'Magneticum', 'TNG_100'] 
orientations = ['xy', 'yz', 'xz']
bands = ['r']

npix_bin = 2 #number of pixels for the binning
factor = 1e12 # scale factor for multiplying the data values
factor_str = format(factor, '.0e')


for sim in sims:
    
    # =============================================================================
    # Paths 
    # =============================================================================
    data_dir = '/Users/z3530031/unsw/Work/Projects/LSST/Chall3-ICL/data/'+sim+'/'
    res_dir = '/Users/z3530031/unsw/Work/Projects/LSST/Chall3-ICL/results/'
    
    #image_file = '00761_0000048_0.05_xy_r_4Mpc'
    
    image_files = glob.glob(data_dir+'*_r_4Mpc.fits') # Load a list 
    
    for clust_files in image_files:
        
        cluster = clust_files[-23:-20]
        orientation = clust_files[-14:-12]
        
        # =============================================================================
        # Create directory to save results
        # =============================================================================
        try:
            # Create target Directory
            os.mkdir(res_dir+sim)
        except:
            pass
        
        try:
            # Create target Directory
            os.mkdir(res_dir+sim+'/'+cluster+'_'+orientation)
        except:
            pass
 
        try:
            # Create target Directory
            os.mkdir(res_dir+sim+'/'+cluster+'_'+orientation+'/Images')
        except:
            pass
        
        
        # =============================================================================
        # Data treatment
        # =============================================================================
        
        for ind, band in enumerate(bands):
        
            # Load images
            im, hdr = fits.getdata(clust_files, header = True)
            wcs = WCS(hdr)
            
            
              
            # =============================================================================
            # Binning data conserving the total flux (HUGE images!)
            # =============================================================================
            bin_im = block_reduce(im, npix_bin) # [sum] provides block summation (and conserves the data sum)
            
                    
            # =============================================================================
            # Cut image to a 2 Mpc radius   
            # =============================================================================
            centre = (int(bin_im.shape[0]/2), int(bin_im.shape[0]/2))
            diam4_mpc = int(bin_im.shape[0]/2)
            
            if diam4_mpc%2==0: # odd num pix. 2Mpc radius
                size_r2mpc = (diam4_mpc+1, diam4_mpc+1)    
            if diam4_mpc%2==1:
                size_r2mpc = (diam4_mpc, diam4_mpc)     
            
            cutout_2Mpc = Cutout2D(bin_im, centre, size_r2mpc, wcs=wcs)
            
            # Check cutout box size over original image
            #plt.imshow(im*factor, origin='lower')
            #cutout_2Mpc.plot_on_original(color='white')
            #plt.show()
        
            
            # =============================================================================
            # Multiply the counts values by a factor 
            # =============================================================================
            new_im = cutout_2Mpc.data*factor
            new_im = np.array(new_im, dtype=float)
            new_hdr = cutout_2Mpc.wcs.to_header()
            
            new_filename = sim+'_'+cluster+'_'+orientation+'_'+band+'_2Mpc'+'_bin'+str(npix_bin)+'_x'+factor_str+'.fits'
            
            
            # Save new .fits file
            save_fits_image(new_im, new_hdr, res_dir+sim+'/'+cluster+'_'+orientation+'/Images/', new_filename)
            
        
            
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 10:52:30 2022

@author: C. Martinez-Lombilla

Script to obtain the BCG+IGL fraction

** Variables are called IGL_... because is easier if any change is required as this 
code is the same as for 5_IGL_Fraction.py but loading diferent images and saving
files in a diferent folder **

BE AWARE of loading the right images.fits

"""

import numpy as np
import os
import glob
import pickle

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.close('all')
plt.ioff()

plt.rcParams["font.family"] = "serif"
plt.rc('text', usetex=False)

from astropy.io import fits
from astropy.table import Table


from CML import LogNorm
from CML import save_fits_image
from CML import sb2counts
from CML import flux_in_ellipse
from CML import luminosity
from CML import dist_ellipse_map
from CML import findRadius
from CML import Pix2Kpc



# =============================================================================
# Initial params
# =============================================================================
sims = ['Horizon_AGN', 'Hydrangea', 'Magneticum', 'TNG_100'] 
clust_ids = {'Horizon_AGN':['048', '078', '157'], 
             'Hydrangea':['047', '087', '116', '213', '294', '296', '299', '338', '353', '727' ], 
             'Magneticum':['002', '033', '114', '140'], 
             'TNG_100':['003', '006', '010']}

orientations = ['xy', 'yz', 'xz']
band = 'r'
pixscale = 0.4

# =============================================================================
# Paths 
# =============================================================================
results_path = '/Users/z3530031/unsw/Work/Projects/LSST/Chall3-ICL/results/'
mask_path = results_path+'Masks/'

table_r = []
table_names = []

for sim in sims:   
    clusters = clust_ids[sim]
    
    for cluster in clusters:
        for orientation in orientations:


            # =============================================================================
            # Path & image
            # =============================================================================
            icl_path = results_path+sim+'/'+cluster+'_'+orientation+'/ICL_fraction/'
            
            image = fits.getdata(icl_path+'Im_IGLfrac-SBlim-26.fits')
            
            
            # =============================================================================
            # Obtain R at wich ICL dominates            
            # =============================================================================
            r_trans_pix = findRadius(image, center = None)
            
            # 1 kpc per arcsec; pixel = 0.4 arcsec
            r_trans_kpc = r_trans_pix*pixscale
            
            # =============================================================================
            #  Save a table with all the info  
            # =============================================================================
            name = sim+'_'+cluster+'_'+orientation
            
            table_r.append(r_trans_kpc)
            table_names.append(name)

table_vec = [table_names, table_r]                    
table_BCG = Table(rows=table_vec, meta={'R_which_ICL_dominates': 'table info'})
table_BCG.write(results_path+'table_R_trans_ICL.fits', overwrite=True)

        

          




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



# =============================================================================
# Initial params
# =============================================================================
sims = ['Horizon_AGN']#, 'Hydrangea', 'Magneticum', 'TNG_100'] 
clust_ids = {'Horizon_AGN':['048']}#, '078', '157'], 'Hydrangea':['213'], 'Magneticum':['114'], 'TNG_100':['006']}
orientations = ['xy']#, 'yz', 'xz']
band = 'r'

pixscale = 0.4
zp = 28.495
sky = 0.00768
sb_lim = 30.0
factor = 1e12

Xc, Yc = [2570,2566] # [5141,5133] H-AGN_48_xy
ell = {'BCG':0}
pa = {'BCG':0}
maxsma = {'BCG':2500}


# =============================================================================
# Paths 
# =============================================================================
results_path = '/Users/z3530031/unsw/Work/Projects/LSST/Chall3-ICL/results/'
mask_path = results_path+'Masks/'


for sim in sims:   
    clusters = clust_ids[sim]
    
    for cluster in clusters:
        for orientation in orientations:
            
            # =============================================================================
            # Prepare diccionaries
            # =============================================================================
            
            TotalFlux = {'Name':'TotalFlux_Lum', 'r':0}
            TotalFlux_err = {'Name':'TotalFluxError', 'r':0}
            
            IGL_2Dmod_sum = {'Name':'BCG+IGL_2Dmod_sum_Lum', 'r':0}
            IGL_2Dmod_sum_err = {'Name':'BCG+IGL_2Dmod_sum_err','r':0}
            IGLfrac_2DMod = {'Name':'BCG+IGLfrac_2DMod','r':0}
            IGLfrac_2DMod_err = {'Name':'BCG+IGLfrac_2DMod_err','r':0}
            
            
            counts_imlim = {'r':0}
            counts_IGLlim = {'r':[]}
            sb_IGLlim = {'r':[]}


            # =============================================================================
            # Create a folder to save script products
            # =============================================================================
            try:
                # Create target Directory
                os.mkdir(results_path+sim+'/'+cluster+'_'+orientation+'/BCG-ICL_fraction')
            except:
                pass
            
            
            # =============================================================================
            # Path BCG+ICL
            # =============================================================================
            icl_path = results_path+sim+'/'+cluster+'_'+orientation+'/BCG-ICL_fraction/'
            

            # =============================================================================
            # 0 - Total masss of the system
            # =============================================================================
            # Load data for total lum
            im_path = results_path+sim+'/'+cluster+'_'+orientation+'/Images/'
            
            im_file = glob.glob(im_path+'*.fits') 
            data = fits.getdata(im_file[0])-sky
            hdr = fits.getheader(im_file[0])

            # SB_lim of the images to counts
            counts_imlim[band] = sb2counts(sb_lim, zp, pixscale)
            
  
            # pdf file to save the plots
            pdf = PdfPages(icl_path+'Apperture_1Mpc.pdf')

            
            # =============================================================================
            # Total flux of the group in counts
            # =============================================================================
            # Mask data below the sb_lim of the image
            data_down = np.copy(data)
            data_down[data<counts_imlim[band]] = np.nan    

            # Get appertures and sum all pixels inside the ellipse 
            Core_sum, Core_sum_err, Core_ap = flux_in_ellipse(data_down, 
                                                          Xc, 
                                                          Yc, 
                                                          maxsma['BCG'], 
                                                          ell['BCG'], 
                                                          pa['BCG'])
        
        
            # =============================================================================
            # Total Flux of the group in counts & luminosity
            # =============================================================================
                
            TotalFlux[band] = Core_sum/factor  
            TotalFlux_err[band] = Core_sum_err/factor  
            #print('Total Flux of the group in counts ('+band+'-band): '+str(TotalFlux))            

        
        
        
            #==============================================================================
            # Plots with the images and apertures
            #==============================================================================
            fig = plt.figure(figsize=(8,8))
            plt.title(sim+'-'+cluster+'-'+orientation, size=18)
            
            # Unmasked data
            plt.imshow(data, origin='lower', cmap='viridis',norm=LogNorm(vmin=0, vmax=5))
            Core_ap.plot(color='white', ls=':')
            
            pdf.savefig(fig)






            # =============================================================================
            # 1 - BCG+IGL fraction
            # =============================================================================
            # Read data
            im_path = results_path+sim+'/'+cluster+'_'+orientation+'/Masks/'
            im_file = glob.glob(im_path+'final*.fits') 
            data = fits.getdata(im_file[0])-sky
            hdr = fits.getheader(im_file[0])        
            
        
            # =============================================================================
            # Total flux within the IGL cicle 
            # =============================================================================
            # Mask data below the sb_lim of the image
            data_down = np.copy(data)
            data_down[data<counts_imlim[band]] = np.nan    
                
            # Get appertures and sum all pixels inside the ellipse 
            IGL_2Dmod_sum[band], IGL_2Dmod_sum_err[band], IGL_ap = flux_in_ellipse(data_down, 
                                                          Xc, 
                                                          Yc, 
                                                          maxsma['BCG'], 
                                                          ell['BCG'], 
                                                          pa['BCG'])
            IGL_2Dmod_sum[band] = IGL_2Dmod_sum[band]/factor  
            IGL_2Dmod_sum_err[band] = IGL_2Dmod_sum_err[band]/factor  
            #print('Total Flux of the IGL in counts ('+band+'-band): '+str(IGL_sum_2Dmod[band]))
            
            # =============================================================================
            # Fraction of IGL - 2D Models & luminosity
            # =============================================================================
            
            IGLfrac_2DMod[band] = IGL_2Dmod_sum[band]/TotalFlux[band]
            IGLfrac_2DMod_err[band] = IGLfrac_2DMod[band] * np.sqrt(np.square(IGL_2Dmod_sum_err[band] / IGL_2Dmod_sum[band]) 
                                                                    + np.square(TotalFlux_err[band] / TotalFlux[band]))
                    
        
            print('***************************************************************')
            print('BCG+IGL fraction from 2D Models ('+band+'-band): '+str(IGLfrac_2DMod[band]))
            print('***************************************************************')
        
            
            # =============================================================================
            #  Save a table with all the info  
            # =============================================================================
            
            table_vec = [TotalFlux, TotalFlux_err, IGL_2Dmod_sum, IGL_2Dmod_sum_err, IGLfrac_2DMod, IGLfrac_2DMod_err]
                    
            table_BCG = Table(rows=table_vec, meta={'BCG+IGL amount and fraction': 'table info'})
            table_BCG.write(icl_path+'BCG-IGLfrac_table.fits', overwrite=True)

    
        pdf.close()
        plt.close('all')
        

          




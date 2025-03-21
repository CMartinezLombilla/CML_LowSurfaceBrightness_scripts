#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 10:52:30 2022

@author: C. Martinez-Lombilla
"""

import numpy as np
import os
import pickle
import glob

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
sims = ['Horizon_AGN', 'Hydrangea', 'Magneticum', 'TNG_100'] 
clust_ids = {'Horizon_AGN':['048', '078', '157'], 
             'Hydrangea':['213', '047', '087', '116', '294', '296', '299', '338', '353', '727' ], 
             'Magneticum':['002', '033', '114', '140'], 
             'TNG_100':['003', '006', '010']}
orientations = ['xy', 'yz', 'xz']
band = 'r'

pixscale = 0.4
zp = 0 #28.495
sky = 0 #0.00768
#sb_lim = 30.3
factor = 1e12


#Xc, Yc = [2570,2566] # [5141,5133] H-AGN_48_xy
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

            IGL_SB_sum = {'Name':'IGL_SB_sum_Lum', 'r':0}
            IGL_SB_sum_err = {'Name':'IGL_SB_sum_err', 'r':0}
            IGLfrac_SB  = {'Name':'IGLfrac_SB', 'r':0}
            IGLfrac_SB_err  = {'Name':'IGLfrac_SB', 'r':0}

            counts_imlim = {'r':0}
            counts_IGLlim = {'r':[]}
            sb_IGLlim = {'r':[]}

  
            # =============================================================================
            # Create a folder to save script products
            # =============================================================================
            try:
                # Create target Directory
                os.mkdir(results_path+sim+'/'+cluster+'_'+orientation+'/ICL_fraction')
            except:
                pass
            
            
            # =============================================================================
            # Path ICL images
            # =============================================================================
            icl_path = results_path+sim+'/'+cluster+'_'+orientation+'/ICL_fraction/'
            
            
            # =============================================================================
            # 0 - Total masss of the system
            # =============================================================================
            # Load data for total lum
            im_path = results_path+sim+'/'+cluster+'_'+orientation+'/Images/'
            
            im_file = glob.glob(im_path+'*.fits') 
            data = (fits.getdata(im_file[0])-sky) /factor
            hdr = fits.getheader(im_file[0])

            # SB_lim of the images to counts
            #counts_imlim[band] = sb2counts(sb_lim, zp, pixscale)
        
            Xc, Yc = np.size(data,1)/2, np.size(data,1)/2 
                
            
            # =============================================================================
            # Total flux of the group in counts
            # =============================================================================
            # Mask data below the sb_lim of the image
            data_down = np.copy(data)
            #data_down[data<counts_imlim[band]] = np.nan    

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
                
            TotalFlux[band] = Core_sum #/factor  
            TotalFlux_err[band] = Core_sum_err #/factor  
            #print('Total Flux of the group in counts ('+band+'-band): '+str(TotalFlux))            

        
        
            ''''
            #==============================================================================
            # Plots with the images and apertures
            #==============================================================================
            # pdf file to save the plots
            pdf = PdfPages(icl_path+'Apperture_1Mpc.pdf')

            fig = plt.figure(figsize=(8,8))
            plt.title(sim+'-'+cluster+'-'+orientation, size=18)
            
            # Unmasked data
            plt.imshow(data, origin='lower', cmap='viridis',norm=LogNorm(vmin=0, vmax=5))
            Core_ap.plot(color='white', ls=':')
            
            pdf.savefig(fig)
            '''


                    
            
            
            
            
            # =============================================================================
            # 1 - IGL fraction from sb limits
            # =============================================================================
            
            sb_IGLlim_list = [26]

            
            # =============================================================================
            # METHOD : convert the sb_lims values to counts
            # =============================================================================
            IGL_sum_SBi = []
            IGL_sum_err_SBi = []
            IGLfrac_SBi = []
            IGLfrac_err_SBi = []
               
                
            for ind, sblim in enumerate(sb_IGLlim_list): # correcting by sbdimming & k-corr
                        
                # IGL sb_lims to counts
                counts_IGLlim[band].append(sb2counts(sblim, zp, pixscale))
                
                # Mask above the sb_limits of the IGL
                data_IGL_ii = np.copy(data_down)
                data_IGL_ii[data>counts_IGLlim[band][ind]] = np.nan
                
                # Get appertures and sum all pixels inside the ellipse 
                IGL_core, IGL_core_err, IGL_ap = flux_in_ellipse(data_IGL_ii, 
                                                              Xc, 
                                                              Yc, 
                                                              maxsma['BCG'], 
                                                              ell['BCG'], 
                                                              pa['BCG'])
                
                
                # =============================================================================
                #  Get appertures and sum all pixels inside the apperture of the core region          
                # =============================================================================                
                IGL_sumii = IGL_core #/factor
                IGL_sum_SBi.append(IGL_sumii)
        
                IGL_errii = IGL_core_err #/factor            
                IGL_sum_err_SBi.append(IGL_errii)
                        
                # =============================================================================
                # Fraction of IGL - at each sb limit
                # =============================================================================     
                
                frac = IGL_sumii/TotalFlux[band]
                IGLfrac_SBi.append(frac)
                
                frac_err = np.sqrt(np.square(IGL_errii / IGL_sumii) 
                                          + np.square(TotalFlux_err[band] / TotalFlux[band]))
                IGLfrac_err_SBi.append(frac_err)
        
                # Save fits image with the fraction of light measured
                save_fits_image(data_IGL_ii, hdr, icl_path, 'Im_IGLfrac-SBlim-'+str(sblim)+'.fits')
        
        
            # guardamos los valores finales en cada banda
            IGL_SB_sum[band] = IGL_sum_SBi
            IGL_SB_sum_err[band] = IGL_sum_err_SBi
            IGLfrac_SB[band] = IGLfrac_SBi  
            IGLfrac_SB_err[band] = IGLfrac_err_SBi          
             
            
            print('***************************************************************')
            print('IGL fraction from SB limits: '+str(IGLfrac_SB))
            print('***************************************************************')
            
            
    
              
    
            
            # =============================================================================
            #  Save a table with all the info  
            # =============================================================================
            
            table_vec = [TotalFlux, TotalFlux_err]
            
            for ii, sblim in enumerate(sb_IGLlim_list):
                sum_sb = {'Name':'IGL_SB_sum '+str(sblim), 'r':IGL_SB_sum['r'][ii]}
                table_vec.append(sum_sb)
                sum_sb_err = {'Name':'IGL_SB_sum_err '+str(sblim),'r':IGL_SB_sum_err['r'][ii]}
                table_vec.append(sum_sb_err)
                
                frac_sb = {'Name':'IGLfrac_SB '+str(sblim), 'r':IGLfrac_SB['r'][ii]}
                table_vec.append(frac_sb)
                frac_sb_err = {'Name':'IGLfrac_SB_err '+str(sblim), 'r':IGLfrac_SB_err['r'][ii]}
                table_vec.append(frac_sb_err)
                
                    
            table_IGL = Table(rows=table_vec, meta={'IGL amount and fraction': 'table info'})
            table_IGL.write(icl_path+'IGLfrac_sbcut_table_NOSKY_NOsblim.fits', overwrite=True)
            
            
            #pdf.close()
            plt.close('all')
            


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
from CML import com
from CML import flux_in_ellipse
from CML import luminosity



# =============================================================================
# Parametros iniciales del grupo/imagen/IGL/aperturas
# =============================================================================
groupID = '400138'
members_IGL = ['BCG', '1660730', '1660646', '1660615']
members_core = ['BCG', '1660730', '1660646', '1660615']
members = ['BCG','1660545', '1660615', '1660646', '1660730', '2301069', 'IGL']

bands = ['g', 'r', 'i']


# Ellipse params for: total flux (Core and IGL), galaxies '1660545' and '2301069'
ell = {'Core':0.52, '1660545':0, '2301069':0}
pa = {'Core':3.33, '1660545':0,'2301069':0}
maxsma = {'Core':500, '1660545':45, '2301069':40}

# =============================================================================
# Prepare diccionaries
# =============================================================================

TotalFlux = {'Name':'TotalFlux', 'g':0, 'r':0, 'i':0}
TotalFlux_err = {'Name':'TotalFluxError', 'g':0, 'r':0, 'i':0}
TotalLum = {'Name':'TotalLum [Lsun]', 'g':0, 'r':0, 'i':0}

IGL_2Dmod_sum = {'Name':'BCG+IGL_2Dmod_sum','g':0, 'r':0, 'i':0}
IGL_2Dmod_sum_err = {'Name':'BCG+IGL_2Dmod_sum_err','g':0, 'r':0, 'i':0}
IGL_2Dmod_Lum = {'Name':'BCG+IGL_2Dmod_Lum [Lsun]', 'g':0, 'r':0, 'i':0}
IGLfrac_2DMod = {'Name':'BCG+IGLfrac_2DMod','g':0, 'r':0, 'i':0}
IGLfrac_2DMod_err = {'Name':'BCG+IGLfrac_2DMod_err','g':0, 'r':0, 'i':0}


counts_imlim = {'g':0, 'r':0, 'i':0}
counts_IGLlim = {'g':[], 'r':[], 'i':[]}
sb_IGLlim = {'g':[], 'r':[], 'i':[]}


# =============================================================================
# Paths
# =============================================================================
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'


# =============================================================================
# Creamos un directorio para guardar los distintos productos del programa
# =============================================================================

try:
    # Create target Directory
    os.mkdir(results_path+'IGL/BCG-IGLfraction')
except:
    pass


# =============================================================================
# Fichero de parametros
# =============================================================================
# Copiamos el fichero de parametros para asegurarnos de que es el mismo

cmd = 'cp '+results_path+'params.py .'

os.system(cmd)
import params

# Parametros del grupo
bcgID = params.bcgID
z = params.z

# Parametros del instrumento y la imagen
zp = params.zp
pixscale = params.pixelscale
SkyRMS = params.SkyRMS
HSC_sb_lim = params.sb_limit_1sigma

# A y Correccion K
k_corrs = params.k_corrs
extinctions = params.extinctions
SkyRMS = params.SkyRMS

# sb dimming: -2.5*np.log10(intens) - 10*np.log10(1+z) -- pq el dimming va como (1+z)**-4 -- Calvi+2014  doi: 10.1088/0004-637X/796/2/102
dimming = 10*np.log10(1+z) 

#Posiciones centrales de cada miembro del grupo, de ls IGL y del core
Xcoords = []
Ycoords = []
for member in members:
    Xcoords.append(int(round(params.Xc[member],0)))
    Ycoords.append(int(round(params.Yc[member],0)))

Xigl = []
Yigl = []
for memberigl in members_IGL:
    Xigl.append(int(round(params.Xc[memberigl],0)))
    Yigl.append(int(round(params.Yc[memberigl],0)))

 
Xcore = []
Ycore = []
for memberc in members_core:
    Xcore.append(int(round(params.Xc[memberc],0)))
    Ycore.append(int(round(params.Yc[memberc],0)))

    
os.system('rm params.py')



# =============================================================================
# Load the 3 masks of the HSC detector
# =============================================================================

HSC_g = fits.getdata(results_path+'masks/intermediate/detectorHSC-g.fits')
HSC_r = fits.getdata(results_path+'masks/intermediate/detectorHSC-r.fits')
HSC_i = fits.getdata(results_path+'masks/intermediate/detectorHSC-i.fits')

HSC_mask = HSC_g + HSC_r + HSC_i






# =============================================================================
# 0 - Total masss of the system
# =============================================================================

# pdf file to save the plots
pdf = PdfPages(results_path+'IGL/BCG-IGLfraction/dist_All_appertures.pdf')

# Find the center of 'light' of the core region
Gmasked = fits.getdata(results_path+'masks/GroupMasked-i.fits')
Xc_core, Yc_core = com(Gmasked, Xcore, Ycore)


for ind, band in enumerate(bands):
    
    # SB_lim of the images to counts
    counts_imlim[band] = sb2counts(HSC_sb_lim[band], zp, pixscale)

 
    # =============================================================================
    # Read data masked, add HSC detector mask, and creat a 'mask' array
    # =============================================================================
    data = fits.getdata(results_path+'masks/GroupMasked-'+band+'.fits')
    
    # Nos aseguramos que estan enmascarados con las 3 mascaras del detector
    data[HSC_mask>0] = np.nan  
    
    mask = np.zeros_like(data)
    mask[np.isnan(data)] = np.nan
    
    # =============================================================================
    # Total flux ellipse (Core)
    # =============================================================================
    # Mask data below the sb_lim of the image
    data_down = np.copy(data)
    data_down[data<counts_imlim[band]] = np.nan    
        
    # Get appertures and sum all pixels inside the ellipse 
    Core_sum, Core_sum_err, Core_ap = flux_in_ellipse(data_down, 
                                                  Xc_core, 
                                                  Yc_core, 
                                                  maxsma['Core'], 
                                                  ell['Core'], 
                                                  pa['Core'],
                                                  SkyRMS[band])
    
        
    
    # =============================================================================
    # Errors: error within the region due to the sky rms       
    # =============================================================================        
    #error_bin=self.rms_sky/np.sqrt(float(np.size(sec_data)))
    
    
    # =============================================================================
    # Member 1660545
    # =============================================================================
    mem_in = members.index('1660545')
    
    # Get appertures and sum all pixels inside the ellipse 
    Mem545_sum, Mem545_sum_err, Mem545_ap = flux_in_ellipse(data_down, 
                                                        Xcoords[mem_in], 
                                                        Ycoords[mem_in], 
                                                        maxsma['1660545'], 
                                                        ell['1660545'], 
                                                        pa['1660545'], 
                                                        SkyRMS[band])
        
    # =============================================================================
    # Member 2301069
    # =============================================================================
    mem_in = members.index('2301069')
    
    # Get appertures and sum all pixels inside the ellipse 
    Mem069_sum, Mem069_sum_err, Mem069_ap = flux_in_ellipse(data_down, 
                                                        Xcoords[mem_in], 
                                                        Ycoords[mem_in], 
                                                        maxsma['2301069'], 
                                                        ell['2301069'], 
                                                        pa['2301069'], 
                                                        SkyRMS[band])
   
    
    # =============================================================================
    # Total Flux of the group in counts & luminosity
    # =============================================================================
        
    TotalFlux[band] = Core_sum + Mem545_sum + Mem069_sum
    TotalFlux_err[band] = np.sqrt(np.square(Core_sum_err) + np.square(Mem545_sum_err) + np.square(Mem069_sum_err))
    #print('Total Flux of the group in counts ('+band+'-band): '+str(TotalFlux))
    TotalLum[band] = luminosity(TotalFlux[band], TotalFlux_err[band], zp, z=z, k_corr=k_corrs[band], band=band)
    

    #==============================================================================
    # Plots with the images and apertures
    #==============================================================================
    fig = plt.figure(ind,figsize=(8,8))
    plt.title(band+'-band; Group ID: '+groupID,size=18)
    
    # Unmasked data
    data_nomask = fits.getdata(results_path+'PSF/PSFcorr-'+band+'-SkySubDist.fits')
    plt.imshow(data_nomask, origin='lower', cmap='viridis',norm=LogNorm(vmin=0, vmax=5))
    Core_ap.plot(color='white', ls=':')
    Mem545_ap.plot(color='white', ls=':')
    Mem069_ap.plot(color='white', ls=':')
    
    # Shadows in the masked regions
    from matplotlib import colors
    cmap = colors.ListedColormap(['k'])
    
    shadow = np.nan_to_num(mask,nan=1)
    shadow[shadow==0]=np.nan
    
    plt.imshow(shadow,cmap=cmap, interpolation='nearest', origin='lower',alpha=0.5)
    plt.xlim(np.shape(data)[0]/2-np.shape(data)[0]/3, np.shape(data)[0]/2+np.shape(data)[0]/3)
    plt.ylim(np.shape(data)[0]/2-np.shape(data)[0]/3, np.shape(data)[0]/2+np.shape(data)[0]/3)
    
    pdf.savefig(fig)

pdf.close()
plt.close('all')
   

    






# =============================================================================
# 1 - BCG+IGL fraction from 2D models [ONLY THIS METHOD]
# =============================================================================

# pdf file to save the plots
pdf = PdfPages(results_path+'IGL/BCG-IGLfraction/dist_BCG-IGL-2DMod_ap'+str(maxsma['Core'])+'pix.pdf')

# Find the center of masss of the BCG+IGL region
Xc_IGL, Yc_IGL = com(Gmasked, Xigl, Yigl)


for ind, band in enumerate(bands):
    
    # =============================================================================
    # Read data and create a mask
    # =============================================================================
    data_IGL = fits.getdata(results_path+'masks/BCG-IGLMasked-'+band+'.fits')
   
    # Nos aseguramos que estan enmascarados con las 3 mascaras del detector
    data_IGL[HSC_mask>0] = np.nan  
    
    mask = np.zeros_like(data_IGL)
    mask[np.isnan(data_IGL)] = np.nan
    
    # =============================================================================
    # Total flux within the IGL ellipse [same as core!]
    # =============================================================================
    
    # Mask data below the sb_lim of the image
    data_down = np.copy(data_IGL)
    data_down[data_IGL<counts_imlim[band]] = np.nan    
        
    # Get appertures and sum all pixels inside the ellipse 
    IGL_2Dmod_sum[band], IGL_2Dmod_sum_err[band], IGL_ap = flux_in_ellipse(data_down, 
                                                  Xc_IGL, 
                                                  Yc_IGL, 
                                                  maxsma['Core'], 
                                                  ell['Core'], 
                                                  pa['Core'],
                                                  SkyRMS[band])
    #print('Total Flux of the IGL in counts ('+band+'-band): '+str(IGL_sum_2Dmod[band]))
    
    # =============================================================================
    # Fraction of IGL - 2D Models & luminosity
    # =============================================================================
    
    IGLfrac_2DMod[band] = IGL_2Dmod_sum[band]/TotalFlux[band]
    IGLfrac_2DMod_err[band] = IGLfrac_2DMod[band] * np.sqrt(np.square(IGL_2Dmod_sum_err[band] / IGL_2Dmod_sum[band]) 
                                                            + np.square(TotalFlux_err[band] / TotalFlux[band]))
    
    IGL_2Dmod_Lum[band] = luminosity(IGL_2Dmod_sum[band], IGL_2Dmod_sum_err[band], zp, z=z, k_corr=k_corrs[band],  band=band)


    print('***************************************************************')
    print('BCG+IGL fraction from 2D Models ('+band+'-band): '+str(IGLfrac_2DMod[band]))
    print('***************************************************************')

    
    # =============================================================================
    #  Plots with the images and apertures
    # =============================================================================

    fig = plt.figure(ind,figsize=(8,8))
    plt.title('Group ID: '+groupID+'; '+band+'-band',size=18)
    
    # Unmasked data
    data_nomask = fits.getdata(results_path+'PSF/PSFcorr-'+band+'-SkySubDist.fits')
    plt.imshow(data_nomask, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=5))
    IGL_ap.plot(color='white', ls=':')
    
    # Shadows in the masked regions
    cmap = colors.ListedColormap(['k'])
    
    shadow = np.nan_to_num(mask,nan=1)
    shadow[shadow==0]=np.nan
    
    plt.imshow(shadow,cmap=cmap, interpolation='nearest', origin='lower',alpha=0.5)
    plt.xlim(np.shape(data_IGL)[0]/2-np.shape(data_IGL)[0]/3, np.shape(data_IGL)[0]/2+np.shape(data_IGL)[0]/3)
    plt.ylim(np.shape(data_IGL)[0]/2-np.shape(data_IGL)[0]/3, np.shape(data_IGL)[0]/2+np.shape(data_IGL)[0]/3)
    
    pdf.savefig(fig)
    
pdf.close()
plt.close('all')


          


# =============================================================================
#  Save a table with all the info  
# =============================================================================

table_vec = [TotalFlux, TotalFlux_err, TotalLum, IGL_2Dmod_sum, IGL_2Dmod_sum_err, IGL_2Dmod_Lum, IGLfrac_2DMod, IGLfrac_2DMod_err]
        
table_BCG = Table(rows=table_vec, meta={'BCG+IGL amount and fraction': 'table info'})
table_BCG.write(results_path+'IGL/BCG-IGLfraction/dist_BCG-IGLfrac_table.fits', overwrite=True)


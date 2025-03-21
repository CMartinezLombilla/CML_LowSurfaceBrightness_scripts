#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 10:52:30 2022

@author: C. Martinez-Lombilla
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



# =============================================================================
# Parametros iniciales del grupo/imagen/IGL/aperturas
# =============================================================================
groupID = '400138'
members_IGL = ['BCG', '1660730', '1660646', '1660615']
members_core = ['BCG', '1660730', '1660646', '1660615']
members = ['BCG','1660545', '1660615', '1660646', '1660730', '2301069', 'IGL']

bands = ['g', 'r', 'i']

# =============================================================================
# Paths
# =============================================================================
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'


# =============================================================================
# Ellipse params for: total flux (Core and IGL), galaxies '1660545' and '2301069'
# =============================================================================

ell_sky = pickle.load(open(results_path+'BCGsky/main_model.p', 'rb' ) )

extra_pa = np.deg2rad(45)
extra_eps = 0.3

ell_core = np.median(ell_sky.isolist.eps)+extra_eps
pa_core = np.median(ell_sky.isolist.pa)+extra_pa

ell = {'Core':ell_core, '1660545':0, '2301069':0}
pa = {'Core':pa_core, '1660545':0,'2301069':0}
maxsma = {'Core':500, '1660545':45, '2301069':40}


# =============================================================================
# Prepare diccionaries
# =============================================================================

TotalFlux = {'Name':'TotalFlux', 'g':0, 'r':0, 'i':0}
TotalFlux_err = {'Name':'TotalFluxError', 'g':0, 'r':0, 'i':0}

IGL_2Dmod_sum = {'Name':'IGL_2Dmod_sum','g':0, 'r':0, 'i':0}
IGL_2Dmod_sum_err = {'Name':'IGL_2Dmod_sum_err','g':0, 'r':0, 'i':0}
IGLfrac_2DMod = {'Name':'IGLfrac_2DMod','g':0, 'r':0, 'i':0}
IGLfrac_2DMod_err = {'Name':'IGLfrac_2DMod_err','g':0, 'r':0, 'i':0}

IGL_SB_sum = {'Name':'IGL_SB_sum','g':0, 'r':0, 'i':0}
IGL_SB_sum_err = {'Name':'IGL_SB_sum_err','g':0, 'r':0, 'i':0}
IGLfrac_SB  = {'Name':'IGLfrac_SB','g':0, 'r':0, 'i':0}
IGLfrac_SB_err  = {'Name':'IGLfrac_SB','g':0, 'r':0, 'i':0}

counts_imlim = {'g':0, 'r':0, 'i':0}
counts_IGLlim = {'g':[], 'r':[], 'i':[]}
sb_IGLlim = {'g':[], 'r':[], 'i':[]}


# =============================================================================
# Creamos un directorio para guardar los distintos productos del programa
# =============================================================================

try:
    # Create target Directory
    os.mkdir(results_path+'IGL/IGLfraction')
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
HSC_sb_lim = params.sb_limit

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
pdf = PdfPages(results_path+'IGL/IGLfraction/All_appertures.pdf')

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
    # Total Flux of the group in counts
    # =============================================================================
        
    TotalFlux[band] = Core_sum + Mem545_sum + Mem069_sum
    TotalFlux_err[band] = np.sqrt(np.square(Core_sum_err) + np.square(Mem545_sum_err) + np.square(Mem069_sum_err))
    #print('Total Flux of the group in counts ('+band+'-band): '+str(TotalFlux))
    
    

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
# 1 - IGL fraction from 2D models
# =============================================================================

# pdf file to save the plots
pdf = PdfPages(results_path+'IGL/IGLfraction/IGL-2DMod_ap'+str(maxsma['Core'])+'pix.pdf')

# Find the center of masss of the IGL region
Xc_IGL, Yc_IGL = com(Gmasked, Xigl, Yigl)


for ind, band in enumerate(bands):
    
    # =============================================================================
    # Read data and create a mask
    # =============================================================================
    data_IGL = fits.getdata(results_path+'masks/IGLMasked-'+band+'.fits')
   
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
    # Fraction of IGL - 2D Models
    # =============================================================================
    
    IGLfrac_2DMod[band] = IGL_2Dmod_sum[band]/TotalFlux[band]
    IGLfrac_2DMod_err[band] = IGLfrac_2DMod[band] * np.sqrt(np.square(IGL_2Dmod_sum_err[band] / IGL_2Dmod_sum[band]) 
                                                            + np.square(TotalFlux_err[band] / TotalFlux[band]))
    
    print('***************************************************************')
    print('IGL fraction from 2D Models ('+band+'-band): '+str(IGLfrac_2DMod[band]))
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
# 2 - IGL fraction from sb limits
# =============================================================================

# UPPER SB LIMITS TO IGL DETECTION (down_sb_lim = HSC_sb_lim)
sb_IGLlim_list = [24, 24.5, 25, 25.5, 26, 26.5, 27, 27.5, 28]

for band in bands:
 
    # =============================================================================
    # Read data, add HSC detector mask, and create 'mask' array
    # =============================================================================
    data_IGL = fits.getdata(results_path+'masks/CoreMasked-'+band+'.fits')
    
    # Add detector masks
    data_IGL[HSC_mask>0] = np.nan  
    
    mask = np.zeros_like(data_IGL)
    mask[np.isnan(data_IGL)] = np.nan
    
    
    # =============================================================================
    # METHOD : convert the sb_lims values to counts
    # =============================================================================
    IGL_sum_SBi = []
    IGL_sum_err_SBi = []
    IGLfrac_SBi = []
    IGLfrac_err_SBi = []
       
    # Masking below image sb_lim (HSC_sb_lim)
    data_down = np.copy(data_IGL)
    data_down[data_IGL<counts_imlim[band]] = np.nan

        
    for ind, sblim in enumerate(sb_IGLlim_list):
                
        # IGL sb_lims to counts
        counts_IGLlim[band].append(sb2counts(sblim, zp, pixscale))
        
        # Mask above the sb_limits of the IGL
        data_IGL_ii = np.copy(data_down)
        data_IGL_ii[data_IGL>counts_IGLlim[band][ind]] = np.nan
        
        # Get appertures and sum all pixels inside the apperture of the core region          
        IGL_sumii, IGL_errii, IGL_ap = flux_in_ellipse(data_IGL_ii, 
                                             Xc_core, 
                                             Yc_core, 
                                             maxsma['Core'], 
                                             ell['Core'], 
                                             pa['Core'],
                                             SkyRMS[band])
        IGL_sum_SBi.append(IGL_sumii)
        IGL_sum_err_SBi.append(IGL_errii)
        
        # =============================================================================
        # Fraction of IGL - at each sb limit
        # =============================================================================     
        frac = IGL_sumii/TotalFlux[band]
        IGLfrac_SBi.append(frac)
        
        frac_err = frac * np.sqrt(np.square(IGL_errii / IGL_sumii) 
                                  + np.square(TotalFlux_err[band] / TotalFlux[band]))
        IGLfrac_err_SBi.append(frac_err)

        # Save fits image with the fraction of light measured
        save_fits_image(data_IGL_ii, 
                        fits.getheader(results_path+'PSF/PSFcorr-'+band+'-SkySubDist.fits'), 
                        results_path+'IGL/IGLfraction/',
                        'IGLfrac-SBlim-'+str(sblim)+'-'+band+'.fits')


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

table_vec = [TotalFlux, TotalFlux_err, IGL_2Dmod_sum, IGL_2Dmod_sum_err, IGLfrac_2DMod, IGLfrac_2DMod_err]

for ii, sblim in enumerate(sb_IGLlim_list):
    sum_sb = {'Name':'IGL_SB_sum '+str(sblim),'g':IGL_SB_sum['g'][ii], 'r':IGL_SB_sum['r'][ii], 'i':IGL_SB_sum['i'][ii]}
    table_vec.append(sum_sb)
    sum_sb_err = {'Name':'IGL_SB_sum_err '+str(sblim),'g':IGL_SB_sum_err['g'][ii], 'r':IGL_SB_sum_err['r'][ii], 'i':IGL_SB_sum_err['i'][ii]}
    table_vec.append(sum_sb_err)

    frac_sb = {'Name':'IGLfrac_SB '+str(sblim),'g':IGLfrac_SB['g'][ii], 'r':IGLfrac_SB['r'][ii], 'i':IGLfrac_SB['i'][ii]}
    table_vec.append(frac_sb)
    frac_sb_err = {'Name':'IGLfrac_SB_err '+str(sblim),'g':IGLfrac_SB_err['g'][ii], 'r':IGLfrac_SB_err['r'][ii], 'i':IGLfrac_SB_err['i'][ii]}
    table_vec.append(frac_sb_err)
    
        
table = Table(rows=table_vec, meta={'IGL amount and fraction': 'table info'})
table.write(results_path+'IGL/IGLfraction/IGLfrac_table.fits', overwrite=True)


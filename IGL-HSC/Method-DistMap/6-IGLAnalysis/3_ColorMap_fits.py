#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 16:49:14 2021

@author: C. Martinez-Lombilla
"""

import os
import numpy as np

import matplotlib.pyplot as plt
plt.ion(),plt.show()
plt.close('all')

from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve

# =============================================================================
# importamos CML
# =============================================================================

import sys
sys.path.append('/Users/felipe/Dropbox/CML/') # path del CML Feli
sys.path.append('/Users/z3530031/Dropbox/CML/') # path del CML Cris

sys.path.append('/Users/felipe/Dropbox/CML/WHT/') # path del fichero de parametros Feli
sys.path.append('/Users/z3530031/Dropbox/CML/WHT/') # path del fichero de parametros Cris



# =============================================================================
# Parametros iniciales del grupo/imagen/IGL/aperturas
# =============================================================================
groupID = '400138'
members = ['IGL']
bands = ['g', 'r', 'i']

# Convolve with a Gaussian kernel?
conv = True
x_sizes = [3, 7, 11, 15]  # 8*sigma+1 is by default
sigma_pix = 3

# IGL-only image or with galaxies?   
image = 'gal' # 'gal' or 'IGL' 



# =============================================================================
# Paths
# =============================================================================
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'

# =============================================================================
# Creamos directorios para guardar los distintos productos del programa
# =============================================================================
try:
    # Create target Directory
    os.mkdir(results_path+'IGL/color/maps')
except:
    pass



# =============================================================================
# Diccionarios para cada uno de los colores
# =============================================================================


colors = ['g-i', 'g-r','r-i']

#files = {'g-i':[results_path+'masks/IGLMasked-g.fits', results_path+'masks/IGLMasked-i.fits'],
         #'g-r':[results_path+'masks/IGLMasked-g.fits', results_path+'masks/IGLMasked-r.fits'], 
         #'r-i':[results_path+'masks/IGLMasked-r.fits', results_path+'masks/IGLMasked-i.fits']}


files = {'gal':{'g-i':[results_path+'PSF/PSFcorr-g-SkySub.fits', results_path+'PSF/PSFcorr-i-SkySub.fits'],
         'g-r':[results_path+'PSF/PSFcorr-g-SkySub.fits', results_path+'PSF/PSFcorr-r-SkySub.fits'], 
         'r-i':[results_path+'PSF/PSFcorr-r-SkySub.fits', results_path+'PSF/PSFcorr-i-SkySub.fits']},
         'IGL':{'g-i':[results_path+'IGLlight-g_unmasked.fits', results_path+'IGLlight-i_unmasked.fits'],
         'g-r':[results_path+'IGLlight-g_unmasked.fits', results_path+'IGLlight-r_unmasked.fits'], 
         'r-i':[results_path+'IGLlight-r_unmasked.fits', results_path+'IGLlight-i_unmasked.fits']}}



for color in colors:    
    # =============================================================================
    # Leemos las imagenes que necesitamos: a - b
    # =============================================================================
    
    band_a = fits.getdata(files[image][color][0])
    band_b = fits.getdata(files[image][color][1])
    
    if conv == True:
        
        # =============================================================================
        # Convolucionamos la mascara con un kernel Gaussiano para agrandarla 
        # =============================================================================
                
        for x_size in x_sizes:
            
            y_size = x_size
            
            kernel= Gaussian2DKernel(sigma_pix, x_size=x_size, y_size=y_size) 
            
            conv_a = convolve(band_a, kernel, preserve_nan=True)
            conv_b = convolve(band_b, kernel, preserve_nan=True)
        
    
            # =============================================================================
            # Mapa de color
            # =============================================================================    
            
            colmap_B = -2.5*np.log10(conv_a)-(-2.5*np.log10(conv_b))
            
            file_name = image+'_Colormap_'+color+'_GF'+str(x_size)


            # =============================================================================
            # Guardamos mapas de color
            # =============================================================================
            header= fits.getheader(files[image][color][0])
            wcs = WCS(header)
            new_hdu = fits.PrimaryHDU(colmap_B)  # we create a PrimaryHDU object to encapsulate the data
            new_hdu.header.update(wcs.to_header())
            new_hdu.writeto(results_path+'IGL/color/maps/'+file_name+'.fits', overwrite=True) # write to a new file
        

 
    if conv == False:
        # =============================================================================
        # Mapa de color
        # =============================================================================    
        
        colmap_B = -2.5*np.log10(band_a)-(-2.5*np.log10(band_b))
        
        file_name = image+'_Colormap_'+color

    
        # =============================================================================
        # Guardamos mapas de color
        # =============================================================================
        header= fits.getheader(files[image][color][0])
        wcs = WCS(header)
        new_hdu = fits.PrimaryHDU(colmap_B)  # we create a PrimaryHDU object to encapsulate the data
        new_hdu.header.update(wcs.to_header())
        new_hdu.writeto(results_path+'IGL/color/maps/'+file_name+'.fits', overwrite=True) # write to a new file
    


import sys
import os
from astropy.visualization import make_lupton_rgb
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel

import numpy as np
import matplotlib.pyplot as plt
plt.ion(),plt.show()
plt.close('all')

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 18

from astropy.table import Table
from astropy.io import fits



# =============================================================================
# Parametros iniciales del grupo/imagen/IGL/aperturas
# =============================================================================
groupID = '400138'
members = ['IGL']
bands = ['g', 'r', 'i']



# =============================================================================
# Paths
# =============================================================================
data_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/data/images/'
press_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Paper/PressReseale-UNSW/'



# =============================================================================
# Diccionarios para cada uno de los archivos
# =============================================================================

files = {'g':data_path+'400138-g.fits',
         'r':data_path+'400138-r.fits', 
         'i':data_path+'400138-i.fits'}


g_data = fits.getdata(files['g'])[950:1650, 850:2100]
r_data = fits.getdata(files['r'])[950:1650, 850:2100]
i_data = fits.getdata(files['i'])[950:1650, 850:2100]


sigma_pix = 1 
kernel= Gaussian2DKernel(sigma_pix)  #X & Y size by default 8*sigma_pix

conv_g = convolve(g_data, kernel)
conv_r = convolve(r_data, kernel)
conv_i = convolve(i_data, kernel)

        
rgb = make_lupton_rgb(i_data, r_data, g_data,
                             Q=11, stretch=0.5, minimum=-0.001,
                             filename=press_path+'rgb_group.pdf')

        
rgb = make_lupton_rgb(conv_i, conv_r, conv_g,
                             Q=13, stretch=0.4, minimum=-0.001,
                             filename=press_path+'rgb_conv_group.pdf')



'''
fig, ax = plt.subplots(figsize=(11,7.5))

plt.imshow(rgb)
        
plt.savefig(press_path+'Group_Color-Lupton.pdf')
'''




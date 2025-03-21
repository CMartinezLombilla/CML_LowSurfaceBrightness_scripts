import os
import numpy as np

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 8})
plt.ion(),plt.show()
plt.close('all')

from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import block_reduce

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



# =============================================================================
# Paths
# =============================================================================
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'
maps_path = results_path+'/IGL/color/maps/'


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


# Parametros del instrumento
zp = params.zp
pixscale = params.pixelscale

# A y Correccion K
k_corrs = params.k_corrs
extinctions = params.extinctions
SkyRMS = params.SkyRMS

# sb dimming: -2.5*np.log10(intens) - 10*np.log10(1+z) -- pq el dimming va como (1+z)**-4 -- Calvi+2014  doi: 10.1088/0004-637X/796/2/102
dimming = 10*np.log10(1+z) 

    
os.system('rm params.py')


# =============================================================================
# Diccionarios para cada uno de los colores
# =============================================================================


colors = ['g-i', 'g-r','r-i']


# BLOCK SIZE IN PIXELS
block = 1

# File names
files = {'g-i':maps_path+'Colormap_B'+str(block)+'_g-i.fits',
         'g-r':maps_path+'Colormap_B'+str(block)+'_g-r.fits', 
         'r-i':maps_path+'Colormap_B'+str(block)+'_r-i.fits'}

# Figure params
color_map = plt.cm.get_cmap('Spectral')
reversed_color_map = color_map.reversed()
fig, axs = plt.subplots(3, 1, sharex=False, sharey=True, figsize=(3,13))
levels = np.array([-0.1, 0.1, 0.3, 0.5, 0.7, 0.9])
row = 0
xmin, xmax = 1000, 1500
ymin, ymax = 1050, 1550


for color in colors:    
  
    # =============================================================================
    # Plots (B2)
    # =============================================================================
    colmap = fits.getdata(files[color])
    Z = colmap

    X = np.arange(colmap.shape[1])
    Y = np.arange(colmap.shape[0])
    
    f = axs[row].contourf(X, Y, Z, levels, cmap=reversed_color_map)
    axs[row].set_xlim([xmin/block,xmax/block])
    axs[row].set_ylim([ymin/block, ymax/block])
    axs[row].set_title(color)

    row+=1

plt.subplots_adjust(left=0.2, right=0.91, hspace=0.5)
fig.colorbar(f, ax=axs, location = 'bottom')
plt.show()
plt.savefig(maps_path+'Colormap-contours-B'+str(block)+'.png', dpi=150)




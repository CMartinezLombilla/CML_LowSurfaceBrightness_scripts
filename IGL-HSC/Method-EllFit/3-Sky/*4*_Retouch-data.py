#*****************************************************************************
#                       RETOUCH MASK IMAGE PROGRAM   [Rix+2004]
#-----------------------------------------------------------------------------
#                             Cristina Martinez 04/21
#*****************************************************************************

# Versions
# %pip install numpy== 1.18.5
# %pip install matplotlib==3.2.2
# %pip install astropy==4.0
# %pip install photutils==0.7.2

#!pip show numpy
#!pip show matplotlib
#!pip show astropy
#!pip show photutils


# Basics
import numpy as np

# Las siguiente dos lineas son para que no pete el MATPLOTLIB fuera de ANACONDA
import matplotlib
matplotlib.use('TkAgg')

# Matplotlib
import matplotlib.pyplot as plt
plt.close('all')
plt.ion()

# Astropy
from astropy.io import fits
from astropy.wcs import WCS


# CML
import sys
sys.path.append('/Users/felipe/Dropbox/CML/')
sys.path.append('/Users/z3530031/Dropbox/CML/')

from CML import HandMask
from CML import Patchs



# =============================================================================
# Cargamos las imagenes y las mascaras 
# =============================================================================

# Group info
groupID = '400138'
band = 'i'  # !!!!!!!! hay que hacer por banda o combinar PSFcorr images ?? 


# Paths & data
masks_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/masks/intermediate/'
load_data = fits.open('/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/PSFcorr-'+band+'.fits')

data = load_data[1].data
data_head = load_data[1].header

#data = fits.getdata(masks_path+'HotColdMasked.fits')
#data_head = fits.getheader(masks_path+'HotColdMasked.fits')
# !!!!!!!! cargar imagen ya enmascarada como hago aqui pero tiene que ser la PSFcorr enmascarada 
 
 
file = masks_path+'HotColdMasked'



#---------------------------------------------------------------
# Primera iteracion
unsh = fits.getdata(masks_path+'UnshapedMasked-retouched.fits')
# Enmacaramos gri con la mascara Hot + Cold + retouched (unsh)
data[np.isnan(unsh)] = np.nan

'''
#---------------------------------------------------------------
# Siguientes iteraciones
retou = fits.getdata(masks_path+'HotColdMasked-retouched.fits')
# Enmacaramos gri con la mascara Hot + Cold + retouched (unsh)
data[np.isnan(retou)] = np.nan
#---------------------------------------------------------------
'''


# =============================================================================
# Retocamos la imagen gri agregando mascaras a mano
# =============================================================================
print('gri image: Add manual masks of a given size over BRIGHT regions around sources')

patchs = Patchs(data,1000)

for i in patchs:
    data[i] = HandMask(data[i], mask_value= np.nan, vmin = 0., vmax=1)



# =============================================================================
# Guardamos la nueva imagen enmascarada 
# =============================================================================

# Create a fits file to save the  mask 
masked_hdu = fits.PrimaryHDU(data)  # we create a PrimaryHDU object to encapsulate the data:
wcs = WCS(data_head)
masked_hdu.header.update(wcs.to_header())
masked_hdu.writeto(file+'-retouched.fits', overwrite=True) # write to a new file


'''
# =============================================================================
# Creamos la mascara final y la guardamos
# =============================================================================

final_mask = np.zeros_like(data)
final_mask[np.isnan(data)] = np.nan

mask_hdu = fits.PrimaryHDU(final_mask) 
wcs = WCS(data_head)
mask_hdu.header.update(wcs.to_header())
mask_hdu.header.append(('HISTORY', ' Final BCG Mask', 'C. Martinez-Lombilla'), end=True)
mask_hdu.writeto('/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/masks/BCGMask.fits', overwrite=True) # write to a new file


masked_hdu = fits.PrimaryHDU(data)  
wcs = WCS(data_head)
masked_hdu.header.update(wcs.to_header())
masked_hdu.writeto('/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/masks/BCGMasked.fits', overwrite=True) # write to a new file

'''














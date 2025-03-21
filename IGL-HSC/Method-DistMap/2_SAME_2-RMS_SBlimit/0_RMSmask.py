#*****************************************************************************
#                       RMS MASK IMAGE PROGRAM - FOR EACH BAND
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

from CML import CircularHandMask
from CML import Patchs



# =============================================================================
# Cargamos las imagenes y las mascaras 
# =============================================================================

# Group info
groupID = '400138'
bands = ['g', 'r', 'i']  

# Path
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'
data_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/data/images/'+groupID





# =============================================================================
# Creamos las mascaras para cada banda sumando la mascara del detector en cada banda 
# =============================================================================


for band in bands:
    
    # PSF mask
    PSFMask = fits.getdata(results_path+'masks/PSFMask-'+band+'.fits')
    #PSFMask = fits.getdata(results_path+'masks/BCGMask-'+band+'.fits')

    HSC_mask = fits.getdata(results_path+'masks/intermediate/detectorHSC-'+band+'.fits')
    data = fits.getdata(results_path+'PSF/PSFcorr-'+band+'.fits')
    #data = fits.getdata(data_path+'-'+band+'.fits')
    
    rms_mask = np.zeros_like(PSFMask)
    rms_mask[np.isnan(PSFMask)] = np.nan
    
    data_masked = data.copy()
    data_masked[np.isnan(rms_mask)] = np.nan
    
    # =============================================================================
    # Retocamos la mascara agregando mascaras a mano si fuera necesario
    # =============================================================================
    print( band+' image: Add manual masks of a given size over BRIGHT regions around sources')
    
    patchs = Patchs(data,2800)
    
    for i in patchs:
        data_masked[i] = CircularHandMask(data_masked[i], mask_value= np.nan, vmin = 0., vmax=0.5)
    
    rms_mask[np.isnan(data_masked)] = np.nan    
    rms_mask[HSC_mask>0] = np.nan    

    # =============================================================================
    # Guardamos la nueva imagen enmascarada 
    # =============================================================================
    
    # Create a fits file to save the  mask 
    masked_hdu = fits.PrimaryHDU(rms_mask)  # we create a PrimaryHDU object to encapsulate the data:
    wcs = WCS(fits.getheader(results_path+'PSF/PSFcorr-'+band+'.fits'))
    masked_hdu.header.update(wcs.to_header())
    masked_hdu.writeto(results_path+'masks/RMSMask-'+band+'.fits', overwrite=True) # write to a new file
    #masked_hdu.writeto(results_path+'masks/BCGMask2-'+band+'.fits', overwrite=True) # write to a new file
    






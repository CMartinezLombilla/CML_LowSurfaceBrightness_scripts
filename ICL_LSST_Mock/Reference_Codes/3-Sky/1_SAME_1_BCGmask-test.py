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
from photutils import detect_threshold
from astropy.convolution import Gaussian2DKernel
from photutils import detect_sources
from photutils import deblend_sources
from astropy.convolution import convolve
from astropy.io import fits
from astropy.wcs import WCS
from regions import read_ds9


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
masks_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/masks/'



# Previous BCG Mask
prevMask = fits.getdata(masks_path+'BCGMask-i.fits')
header = fits.getheader(masks_path+'BCGMask-i.fits')


# Data
gri = fits.getdata('/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/data/images/'+groupID+'-gri.fits')


'''

# Core Mask
coreMask = fits.getdata(masks_path+'CoreMask.fits')
header = fits.getheader(masks_path+'CoreMask.fits')


# =============================================================================
# Read regions to mask (all core galaxies but the BCG) within the bigCore region
# =============================================================================
regCore = read_ds9(masks_path+'ds9Regions/BigCoreMask.reg')

regsCoreMask = np.zeros_like(prevMask)
shape = regsCoreMask.shape
wcs = WCS(header)

for i in regCore:
    regsCoreMask = regsCoreMask +i.to_pixel(wcs).to_mask('center').to_image(shape)



regGlxs = read_ds9(masks_path+'ds9Regions/CoreGroupGal.reg')
regGlxsMask = np.zeros_like(prevMask)

for i in regGlxs:
    regGlxsMask = regGlxsMask +i.to_pixel(wcs).to_mask('center').to_image(shape)
    
    

# =============================================================================
# Creamos las mascaras para cada banda sumando la mascara del core y de las galaxias del grupo excepto la BCG
# =============================================================================


UnmaskCore = fits.getdata(masks_path+'intermediate/unmaskedCore-i.fits')

BCG_mask = np.zeros_like(UnmaskCore)
BCG_mask[np.isnan(UnmaskCore)] = np.nan


BCG_mask[regsCoreMask!=0] = coreMask[regsCoreMask!=0]
BCG_mask[regGlxsMask!=0] = np.nan

'''


# =============================================================================
# Convolved image
# =============================================================================
# Convolve the image with a gaussian kernel 
sigma_pix = 3
kernel= Gaussian2DKernel(sigma_pix) #X & Y size by default 8*sigma_pix
conv_gri = convolve(gri, kernel)


# Masking with the previous BCG mask
prev_mask = np.isnan(prevMask)

data = np.zeros_like(conv_gri)
data[prev_mask] = np.nan
data[~prev_mask] = conv_gri[~prev_mask]



# =============================================================================
# Hand masking
# =============================================================================
patchs = Patchs(data,900) 

handmask = np.zeros_like(conv_gri)


for i in patchs:
	handmask[i] = CircularHandMask(data[i], mask_value=np.nan, vmin=-0.001, vmax=0.2)



# =============================================================================
# Add mask in the core region to the BCGmask in each band with the HSC-detector mask
# =============================================================================

for band in bands:
    
    BCGMask = fits.getdata(masks_path+'BCGMask-'+band+'.fits')
    header = fits.getheader(masks_path+'BCGMask-'+band+'.fits')


    # Add detector mask from the core region
    HSC_mask = fits.getdata(masks_path+'intermediate/detectorHSC-'+band+'.fits')   
    BCGMask[HSC_mask>0] = np.nan  
    
    # Add manual masks we just made
    BCGMask[np.isnan(handmask)] = np.nan  

     
    # =============================================================================
    # Guardamos la nueva imagen mascara de la BCG 
    # =============================================================================
    
    # Create a fits file to save the  mask 
    masked_hdu = fits.PrimaryHDU(BCGMask)  # we create a PrimaryHDU object to encapsulate the data:
    wcs = WCS(header)
    masked_hdu.header.update(wcs.to_header())
    masked_hdu.writeto(masks_path+'BCGMask3-'+band+'.fits', overwrite=True) # write to a new file
    


    # =============================================================================
    # Guardamos la nueva imagen enmascarada 
    # =============================================================================
    data = fits.getdata(results_path+'PSF/PSFcorr-'+band+'.fits')
    data[np.isnan(BCGMask)] = np.nan 
    # Create a fits file to save the  mask 
    masked_hdu = fits.PrimaryHDU(data)  # we create a PrimaryHDU object to encapsulate the data:
    wcs = WCS(header)
    masked_hdu.header.update(wcs.to_header())
    masked_hdu.writeto(masks_path+'BCGMasked3-'+band+'.fits', overwrite=True) # write to a new file
    










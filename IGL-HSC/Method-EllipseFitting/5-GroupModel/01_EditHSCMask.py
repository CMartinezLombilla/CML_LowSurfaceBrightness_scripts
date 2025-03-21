

import matplotlib.pyplot as plt

from astropy.io import fits
from regions import read_ds9
from astropy.wcs import WCS
import numpy as np
from photutils.aperture import RectangularAperture


# =============================================================================
# Info inicial
# =============================================================================

# Group info
groupID = '400138'


# **************************************
# **************************************

# What do I want to do?
Task = 'unmask'                   # mask or unmask?
In_maskFileName = 'CoreMask'    # 'GroupMask' ; 
Out_maskFileName = 'CoreMaskEdit'   # 'CoreMask' ; 'GroupMask'

# Region size and center to unmask
Xc = 1230
Yc = 1282
w_reg = 200
h_reg = 15


# **************************************
# **************************************


# Bands
bands = ['g']#, 'r', 'i']   # !!!!!!!! hay que hacer por banda o combinar PSFcorr images ?? 

# Path
masks_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/masks/'
data_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/PSF/'


        
# =============================================================================
# Rectangular aperture to subtract    
# =============================================================================
region = RectangularAperture([Xc, Yc], w_reg, h_reg)



for band in bands:
    # =============================================================================
    # Cargamos las mascaras 
    # =============================================================================
    
    mask = fits.getdata(masks_path+In_maskFileName+'-'+band+'.fits')
    header = fits.getheader(masks_path+In_maskFileName+'-'+band+'.fits')
    wcs = WCS(header)
    shape = mask.shape
    
        
    # =============================================================================
    # Rectangular aperture to subtract    
    # =============================================================================
    im_region = region.to_mask('center').to_image(shape)


    # =============================================================================
    # Enmascaramos o desenmascaramos al gusto   
    # =============================================================================
  
    if Task == 'mask':
        mask[im_region!=0]=np.nan

    
    if Task == 'unmask':
                
        mask[im_region!=0]=0

    
        
    # =============================================================================
    # Combinamos con la mascara del detector
    # =============================================================================    
    
    new_hdu = fits.PrimaryHDU(mask)  # we create a PrimaryHDU object to encapsulate the data
    new_hdu.header.update(wcs.to_header())
    new_hdu.writeto(masks_path+Out_maskFileName+'-'+band+'.fits', overwrite=True) # write to a new file
    
     
    # Guardamos imagen enmascarada
    masked = fits.getdata(data_path+'PSFcorr-'+band+'-SkySub.fits')    
    masked[np.isnan(mask)]=np.nan  
    
    new_hdu = fits.PrimaryHDU(masked)  # we create a PrimaryHDU object to encapsulate the data
    new_hdu.header.update(wcs.to_header())
    new_hdu.writeto(masks_path+Out_maskFileName+'ed-'+band+'.fits', overwrite=True) # write to a new file
  
 
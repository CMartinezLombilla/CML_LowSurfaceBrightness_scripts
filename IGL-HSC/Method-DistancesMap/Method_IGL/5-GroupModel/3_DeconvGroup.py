

import matplotlib.pyplot as plt

from astropy.io import fits
from regions import read_ds9
from astropy.wcs import WCS
import numpy as np

from CML import save_fits_image


# =============================================================================
# Input info
# =============================================================================

# Group info
groupID = '400138'
bands = ['g','r','i']   # !!!!!!!! hay que hacer por banda o combinar PSFcorr images ?? 
members = ['1660545', '2301069', '1660615']  # Select/coment members as needed 

# Data path
masks_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/masks/'

# DS9 Regs over the galaxy members outside the group core
regMem = []
for mem in members:
    regMem.append(read_ds9(masks_path+'ds9Regions/'+mem+'.reg'))


# =============================================================================
# Actual process
# =============================================================================
    
for band in bands:  

    # Models path, Masks and data
    models_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/2D-Imfit/'+band+'/distMap/'

    Gmask = fits.getdata(masks_path+'GroupMasked-'+band+'.fits')    
    header = fits.getheader(masks_path+'GroupMask-'+band+'.fits')
    
    
    deconv_mod = fits.getdata(models_path+'2DmodSers-'+band+'.fits')
    conv_res = fits.getdata(models_path+'2DresSersPSF-'+band+'.fits')
    
    # =============================================================================
    # Combine images: PSFDeconvMod + ResPSFConvMod
    # =============================================================================
    deconv_group = deconv_mod + conv_res
    
    masked = np.copy(deconv_group)
    masked[np.isnan(Gmask)]=np.nan
    
    # =============================================================================
    # Add data to the masked galaxies outside the core of the group    
    # =============================================================================

    imRegMem = np.zeros_like(Gmask)   
    shape = imRegMem.shape
    wcs = WCS(header)
    
    for mem in range(len(members)):
        for i in regMem[mem]:    
            imRegMem = imRegMem +i.to_pixel(wcs).to_mask('center').to_image(shape)

    masked[imRegMem!=0]=Gmask[imRegMem!=0]
    
    
    # =============================================================================
    # Save masked image
    # =============================================================================
    save_fits_image(masked, header, masks_path, 'DeconvGroupMasked-'+band+'.fits')
    
    
    



import matplotlib.pyplot as plt

from astropy.io import fits
from regions import read_ds9
from astropy.wcs import WCS
import numpy as np


# =============================================================================
# Cargamos las imagenes y las mascaras 
# =============================================================================

# Group info
groupID = '400138'
bands = ['g','r','i']   # !!!!!!!! hay que hacer por banda o combinar PSFcorr images ?? 

# Path
masks_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/masks/'
data_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/PSF/'


for band in bands:  

    # Mascara a desenmascarar el core
    Gmask = fits.getdata(masks_path+'GroupMask-'+band+'.fits')    
    header = fits.getheader(masks_path+'GroupMask-'+band+'.fits')
    
    maskHSC = fits.getdata(masks_path+'intermediate/detectorHSC-'+band+'.fits')

    # =============================================================================
    # Leemos las regiones a desenmascarar
    # =============================================================================
    reg = read_ds9(masks_path+'ds9Regions/TinyThingsCore.reg')
    
    
    # =============================================================================
    # Hacemos un bucle que pasa por las regiones para crear la máscara
    # =============================================================================
    
    mask = np.zeros_like(Gmask)
    shape = mask.shape
    wcs = WCS(header)
    
    for i in reg:    
       mask = mask +i.to_pixel(wcs).to_mask('center').to_image(shape)
       
    
    # =============================================================================
    # Enmascaramos
    # =============================================================================
    Gmask[np.isnan(maskHSC)]=0
    Gmask[mask!=0]=np.nan

    # =============================================================================
    # Guardamos la nueva mascara y la imagen enmascarada
    # =============================================================================
    new_hdu = fits.PrimaryHDU(Gmask)  # we create a PrimaryHDU object to encapsulate the data
    wcs = WCS(header)
    new_hdu.header.update(wcs.to_header())
    new_hdu.writeto(masks_path+'GroupMask-'+band+'.fits', overwrite=True) # write to a new file

    masked = fits.getdata(data_path+'PSFcorr-'+band+'-SkySubDist.fits')    
    masked[Gmask!=0]=np.nan
    new_hdu = fits.PrimaryHDU(masked)  # we create a PrimaryHDU object to encapsulate the data
    wcs = WCS(header)
    new_hdu.header.update(wcs.to_header())
    new_hdu.writeto(masks_path+'GroupMasked-'+band+'.fits', overwrite=True) # write to a new file


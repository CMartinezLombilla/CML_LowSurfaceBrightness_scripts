

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
bands = ['g', 'r', 'i']   # !!!!!!!! hay que hacer por banda o combinar PSFcorr images ?? 

# Path
masks_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/masks/'


# =============================================================================
# Leemos las regiones a desenmascarar
# =============================================================================
reg = read_ds9(masks_path+'ds9Regions/BigCoreMask.reg')
  


for band in bands:
    # =============================================================================
    # Cargamos las mascaras 
    # =============================================================================
    
    # Mascara usada para calcular el RMS de cada banda
    maskRMS = fits.getdata(masks_path+'RMSMask-'+band+'.fits')
    header = fits.getheader(masks_path+'RMSMask-'+band+'.fits')
  
    
    # =============================================================================
    # Hacemos un bucle que pasa por las regiones para crear la máscara
    # =============================================================================
    
    unmask = np.zeros_like(maskRMS)
    shape = unmask.shape
    wcs = WCS(header)
    
    for i in reg:    
       unmask = unmask +i.to_pixel(wcs).to_mask('center').to_image(shape)
       
    
    # =============================================================================
    # Desenmascaramos
    # =============================================================================
    maskRMS[unmask!=0]=0   
    
    
    # =============================================================================
    # Combinamos con la mascara del detector
    # =============================================================================
    maskHSC = fits.getdata(masks_path+'intermediate/detectorHSC-'+band+'.fits')
    maskRMS[maskHSC!=0]=np.nan   
    
    
    new_hdu = fits.PrimaryHDU(maskRMS)  # we create a PrimaryHDU object to encapsulate the data
    wcs = WCS(header)
    new_hdu.header.update(wcs.to_header())
    new_hdu.writeto(masks_path+'intermediate/unmaskedCore-'+band+'.fits', overwrite=True) # write to a new file
    

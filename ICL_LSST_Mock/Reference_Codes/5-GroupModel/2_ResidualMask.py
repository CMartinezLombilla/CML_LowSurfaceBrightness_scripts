

from astropy.io import fits
from regions import read_ds9
from astropy.wcs import WCS
import numpy as np


# =============================================================================
# Cargamos las imagenes y las mascaras 
# =============================================================================

# Group info
groupID = '400138'
bands = ['g', 'r', 'i'] 


# Path
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'
masks_path = results_path+'masks/'


# =============================================================================
# Mascaras previas para aplicar a todas las bandas por igual 
# =============================================================================    
Cmaskg = fits.getdata(masks_path+'CoreMask-g.fits')
Cmaskr = fits.getdata(masks_path+'CoreMask-r.fits')
Cmaski = fits.getdata(masks_path+'CoreMask-i.fits')

Cmask = Cmaskg + Cmaskr + Cmaski


for band in bands:
          
    # Data
    data = fits.getdata(results_path+'IGLlight-'+band+'.fits')
    header = fits.getheader(results_path+'IGLlight-'+band+'.fits')
    
    # =============================================================================
    # Leemos las regiones a enmascarar
    # =============================================================================
    reg = read_ds9(masks_path+'ds9Regions/ResidualsSersicMask.reg')
    
    
    # =============================================================================
    # Hacemos un bucle que pasa por las regiones para crear la máscara
    # =============================================================================
    
    mask = np.zeros_like(Cmask)
    shape = mask.shape
    wcs = WCS(header)
    
    for i in reg:    
       mask = mask +i.to_pixel(wcs).to_mask('center').to_image(shape)
       
    
    # =============================================================================
    # Enmascaramos las regiones y la mascara del detector
    # =============================================================================
    Cmask[mask!=0]=np.nan   
    
    # =============================================================================
    # Guardamos la nueva mascara + imagen enmascarada
    # =============================================================================
    new_hdu = fits.PrimaryHDU(Cmask)  # we create a PrimaryHDU object to encapsulate the data
    wcs = WCS(header)
    new_hdu.header.update(wcs.to_header())
    new_hdu.writeto(masks_path+'IGLMask-'+band+'.fits', overwrite=True) # write to a new file
    
    data[np.isnan(Cmask)] = np.nan  
    new_hdu2 = fits.PrimaryHDU(data)  # we create a PrimaryHDU object to encapsulate the data
    wcs = WCS(header)
    new_hdu2.header.update(wcs.to_header())
    new_hdu2.writeto(masks_path+'Dist_IGLMasked-'+band+'.fits', overwrite=True) # write to a new file
    



from astropy.io import fits
from regions import read_ds9
from astropy.wcs import WCS
import numpy as np

from CML import save_fits_image


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
#Cmaskg = fits.getdata(masks_path+'CoreMask-g.fits')
#Cmaskr = fits.getdata(masks_path+'CoreMask-r.fits')
#Cmaski = fits.getdata(masks_path+'CoreMask-i.fits')

#Cmask = Cmaskg + Cmaskr + Cmaski


for band in bands:
          
    # Data
    igl = fits.getdata(results_path+'pretty_IGLMasked-'+band+'.fits')
    #bcg = fits.getdata(results_path+'BCG_and_IGLlight-'+band+'.fits')
    
    header = fits.getheader(results_path+'pretty_IGLMasked-'+band+'.fits')
    
    # =============================================================================
    # Leemos las regiones a enmascarar
    # =============================================================================
    reg_igl = read_ds9(masks_path+'ds9Regions/ResidualsSersicMask.reg')    
    #reg_bcg = read_ds9(masks_path+'ds9Regions/BCGfracMask.reg')
    
    # =============================================================================
    # Hacemos un bucle que pasa por las regiones para crear la máscara
    # =============================================================================
    
    mask = np.zeros_like(igl)
    shape = mask.shape
    wcs = WCS(header)
    
    for i in reg_igl:    
       mask = mask +i.to_pixel(wcs).to_mask('center').to_image(shape)
    
    # Enmascaramos las regiones y la mascara del detector
    igl_mask = np.copy(mask)
    igl_mask[mask!=0]=np.nan   

    # -------------------------------------------------------------------------
    '''
    mask = np.zeros_like(Cmask)
    
    for i in reg_bcg:    
       mask = mask +i.to_pixel(wcs).to_mask('center').to_image(shape)
    
    # Enmascaramos las regiones y la mascara del detector
    bcg_mask = np.copy(Cmask)
    bcg_mask[mask!=0]=np.nan  
    '''
    # =============================================================================
    # Guardamos la nueva mascara + imagen enmascarada
    # =============================================================================    
    # IGL
    #save_fits_image(igl_mask, header, masks_path, 'IGLMask-'+band+'.fits')

    igl[np.isnan(igl_mask)] = np.nan
    save_fits_image(igl, header, results_path, 'press_IGLMasked-'+band+'.fits')
    '''
    # BCG+IGL
    bcg[np.isnan(bcg_mask)] = np.nan
    save_fits_image(bcg, header, masks_path, 'BCG-IGLMasked-'+band+'.fits')
    ''' 

    

    

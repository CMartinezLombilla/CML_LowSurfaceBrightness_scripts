

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
members = ['1660545', '2301069', '1660615']
membersForeg = ['1660545', '2301069', '1660615']

bands = ['g', 'r', 'i']   # !!!!!!!! hay que hacer por banda o combinar PSFcorr images ?? 

# Path
masks_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/masks/'


# =============================================================================
# Leemos las regiones a desenmascarar
# =============================================================================
regCoreGroup = read_ds9(masks_path+'ds9Regions/CoreGroupGal.reg')

regMem = []
for mem in members:
    regMem.append(read_ds9(masks_path+'ds9Regions/'+mem+'.reg'))

regMemForeg = []
for foreg in membersForeg:
    regMemForeg.append(read_ds9(masks_path+'ds9Regions/'+foreg+'-foreg.reg'))


for band in bands:
    # =============================================================================
    # Cargamos las mascaras 
    # =============================================================================
    
    # Mascara con la BCG desenmascarada de cada banda
    maskBCG = fits.getdata(masks_path+'BCGMask-'+band+'.fits')
    header = fits.getheader(masks_path+'BCGMask-'+band+'.fits')
  
    
    # =============================================================================
    # Desenmascaramos: Hacemos un bucle que pasa por las regiones para crear la máscara
    # =============================================================================
    
    unmask = np.zeros_like(maskBCG)
    shape = unmask.shape
    wcs = WCS(header)
    
    for i in regCoreGroup:    
       unmask = unmask +i.to_pixel(wcs).to_mask('center').to_image(shape)
    
    for mem in range(len(members)):
        for i in regMem[mem]:    
            unmask = unmask +i.to_pixel(wcs).to_mask('center').to_image(shape)
    

    maskBCG[unmask!=0]=0   
    

    # =============================================================================
    # Enmascaramos los foreground objects
    # =============================================================================
 
    mask = np.zeros_like(maskBCG)
    shape = unmask.shape
    wcs = WCS(header)
    
    for memF in range(len(membersForeg)):
        for i in regMemForeg[memF]:    
            mask = mask +i.to_pixel(wcs).to_mask('center').to_image(shape)
    
    maskBCG[mask!=0]=np.nan   
    

    
    # =============================================================================
    # Combinamos con la mascara del detector
    # =============================================================================
    maskHSC = fits.getdata(masks_path+'intermediate/detectorHSC-'+band+'.fits')
    maskBCG[maskHSC!=0]=np.nan   
    
    
    new_hdu = fits.PrimaryHDU(maskBCG)  # we create a PrimaryHDU object to encapsulate the data
    wcs = WCS(header)
    new_hdu.header.update(wcs.to_header())
    new_hdu.writeto(masks_path+'intermediate/GroupMask-'+band+'.fits', overwrite=True) # write to a new file
    

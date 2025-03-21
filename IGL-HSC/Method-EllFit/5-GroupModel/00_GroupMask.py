

import matplotlib.pyplot as plt

from astropy.io import fits
from regions import read_ds9
from astropy.wcs import WCS
import numpy as np


# =============================================================================
# Info inicial
# =============================================================================

# Group info
groupID = '400138'
members = ['1660545', '2301069', '1660615']  # Select/coment members as needed 
membersForeg = ['1660545', '2301069', '1660615']

# **************************************
# **************************************

# What do I want to do?
Task = 'mask'                   # mask or unmask?
Core = False                    # mask of the core galaxies
GroupMem = True                 # mask of each galaxy member out-core 
In_maskFileName = 'GroupMask'   # 'GroupMask' ; 
Out_maskFileName = 'CoreMask'   # 'CoreMask' ; 'GroupMask'

# **************************************
# **************************************


# Bands
bands = ['g', 'r', 'i']   # !!!!!!!! hay que hacer por banda o combinar PSFcorr images ?? 

# Path
masks_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/masks/'
data_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/PSF/'




# =============================================================================
# Leemos las regiones a enmascarar/desenmascarar de cada miembro
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
    # Cargamos las mascaras que enmascaran todo menos los miembros del grupo
    # =============================================================================
    
    maskGroup = fits.getdata(masks_path+In_maskFileName+'-'+band+'.fits')
    header = fits.getheader(masks_path+In_maskFileName+'-'+band+'.fits')
    
    # =============================================================================
    # Hacemos un bucle que pasa por las regiones para crear cada imagen
    # =============================================================================
    
    imRegCore = np.zeros_like(maskGroup)
    imRegMem = np.zeros_like(maskGroup)
    
    shape = imRegCore.shape
    wcs = WCS(header)
    
    for i in regCoreGroup:    
       imRegCore = imRegCore +i.to_pixel(wcs).to_mask('center').to_image(shape)
    
    for mem in range(len(members)):
        for i in regMem[mem]:    
            imRegMem = imRegMem +i.to_pixel(wcs).to_mask('center').to_image(shape)
    

    # =============================================================================
    # Y lo mismo con los foreground objects
    # =============================================================================
 
    imRegMemF = np.zeros_like(maskGroup)
    
    for memF in range(len(membersForeg)):
        for i in regMemForeg[memF]:    
            imRegMemF = imRegMemF +i.to_pixel(wcs).to_mask('center').to_image(shape)
            
            
            
    # =============================================================================
    # Enmascaramos o desenmascaramos al gusto   
    # =============================================================================
  
    if Task == 'mask':
        
        if Core == True:
            # Core
            maskGroup[imRegCore!=0]=np.nan
        
        if GroupMem == True:
            # Group members
            maskGroup[imRegMem!=0]=np.nan

    
    if Task == 'unmask':
                
        if Core == True:
            # Core
            maskGroup[imRegCore!=0]=0
        
        if GroupMem == True:   
            # Group members
            maskGroup[imRegMem!=0]=0
            # Foreg. sources around Group members
            maskGroup[imRegMemF!=0]=np.nan       
    
        
    # =============================================================================
    # Combinamos con la mascara del detector
    # =============================================================================
    maskHSC = fits.getdata(masks_path+'intermediate/detectorHSC-'+band+'.fits')
    maskGroup[maskHSC!=0]=np.nan   
    
    new_hdu = fits.PrimaryHDU(maskGroup)  # we create a PrimaryHDU object to encapsulate the data
    wcs = WCS(header)
    new_hdu.header.update(wcs.to_header())
    new_hdu.writeto(masks_path+Out_maskFileName+'-'+band+'.fits', overwrite=True) # write to a new file
    
    
    # Guardamos imagen enmascarada    
    masked = fits.getdata(data_path+'PSFcorr-'+band+'-SkySub.fits')    
    masked[np.isnan(maskGroup)]=np.nan  
    
    new_hdu = fits.PrimaryHDU(masked)  # we create a PrimaryHDU object to encapsulate the data
    new_hdu.header.update(wcs.to_header())
    new_hdu.writeto(masks_path+Out_maskFileName+'ed-'+band+'.fits', overwrite=True) # write to a new file
  
    
  


import matplotlib.pyplot as plt

from astropy.io import fits
from regions import read_ds9
from astropy.wcs import WCS
import numpy as np

from CML import HandMask
from CML import Patchs

# =============================================================================
# Cargamos las imagenes y las mascaras 
# =============================================================================

# Group info
groupID = '400138'
band = 'i'   # !!!!!!!! hay que hacer por banda o combinar PSFcorr images ?? 
max_gal_numb = 4


# Path
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'
masks_path = results_path+'masks/'

# Data
data = fits.getdata(results_path+'GalaxySubtraction/'+band+'-gal'+str(max_gal_numb)+'.fits')
header = fits.getheader(results_path+'GalaxySubtraction/'+band+'-gal'+str(max_gal_numb)+'.fits')

# Previous Core mask
core_mask = fits.getdata(masks_path+'CoreMask.fits')

# =============================================================================
# Leemos las regiones a enmascarar
# =============================================================================
reg = read_ds9(masks_path+'ds9Regions/CoreGalMasks.reg')


# =============================================================================
# Hacemos un bucle que pasa por las regiones para crear la máscara
# =============================================================================

mask = np.zeros_like(core_mask)
shape = mask.shape
wcs = WCS(header)

for i in reg:    
   mask = mask +i.to_pixel(wcs).to_mask('center').to_image(shape)
   

# =============================================================================
# Desenmascaramos
# =============================================================================
core_mask[mask!=0]=np.nan   


# =============================================================================
# Ponemos mascaras alrededor del core para enmascarar completmante todas las fuentes
# =============================================================================

data[np.isnan(core_mask)] = np.nan

# Retocamos la imagen gri agregando mascaras a mano
print('gri image: Add manual masks of a given size over BRIGHT regions around sources')

patchs = Patchs(data,500)

for i in patchs:
    data[i] = HandMask(data[i], mask_value= np.nan, vmin = 0., vmax=1)



# =============================================================================
# Guardamos la nueva mascara
# =============================================================================
core_mask[np.isnan(data)]=np.nan  

new_hdu = fits.PrimaryHDU(core_mask)  # we create a PrimaryHDU object to encapsulate the data
wcs = WCS(header)
new_hdu.header.update(wcs.to_header())
new_hdu.writeto(results_path+'CoreMask.fits', overwrite=True) # write to a new file


new_hdu2 = fits.PrimaryHDU(data)  # we create a PrimaryHDU object to encapsulate the data
new_hdu2.header.update(wcs.to_header())
new_hdu2.writeto(results_path+'CoreMasked.fits', overwrite=True) # write to a new file



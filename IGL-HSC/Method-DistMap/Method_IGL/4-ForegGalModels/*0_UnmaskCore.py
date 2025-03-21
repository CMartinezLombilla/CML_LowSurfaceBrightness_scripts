

from astropy.io import fits
from regions import read_ds9
from astropy.wcs import WCS
import numpy as np


# =============================================================================
# Cargamos las imagenes y las mascaras 
# =============================================================================

# Group info
groupID = '400138'
band = 'i'   # !!!!!!!! hay que hacer por banda o combinar PSFcorr images ?? 

# Path
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'
masks_path = results_path+'masks/'

# Data
header = fits.getheader(results_path+'PSF/PSFcorr-'+band+'.fits')


# Mascara a desenmascarar el core
BCGmask = fits.getdata(masks_path+'BCGmask.fits')


# =============================================================================
# Leemos las regiones a desenmascarar
# =============================================================================
reg = read_ds9(masks_path+'ds9Regions/CoreUnmask.reg')


# =============================================================================
# Hacemos un bucle que pasa por las regiones para crear la máscara
# =============================================================================

unmask = np.zeros_like(BCGmask)
shape = unmask.shape
wcs = WCS(header)

for i in reg:    
   unmask = unmask +i.to_pixel(wcs).to_mask('center').to_image(shape)
   

# =============================================================================
# Desenmascaramos
# =============================================================================
BCGmask[unmask!=0]=0   

# =============================================================================
# Guardamos la nueva mascara
# =============================================================================
new_hdu = fits.PrimaryHDU(BCGmask)  # we create a PrimaryHDU object to encapsulate the data
wcs = WCS(header)
new_hdu.header.update(wcs.to_header())
new_hdu.writeto(masks_path+'CoreMask.fits', overwrite=True) # write to a new file


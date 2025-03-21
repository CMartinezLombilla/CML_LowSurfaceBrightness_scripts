import matplotlib.pyplot as plt

from astropy.io import fits
from regions import read_ds9
from astropy.wcs import WCS
import numpy as np

# =============================================================================
# Cargamos los datos 
# =============================================================================
mask = fits.getdata('masks/PSFMask.fits')
shape = mask.shape
header = fits.getheader('masks/PSFMask.fits')
wcs = WCS(header)


reg = read_ds9('ds9Regions/Unmask.reg')


# =============================================================================
# Hacemos un bucle para crear la máscara
# =============================================================================

region = np.zeros_like(mask)
for i in reg:    
   region = region +i.to_pixel(wcs).to_mask('center').to_image(shape)
mask[region !=0] = 0
mask = np.nan_to_num(mask,nan=1)
# =============================================================================
#  Sumamos el resultado a la mascaras finales
# =============================================================================

lista = ['masks/mask_s2.fits','masks/mask_s3.fits','masks/mask_s1.fits']
    
for i in lista:
    final_mask = fits.getdata(i) + mask
    final_mask[final_mask!=0]=np.nan
    new_hdu = fits.PrimaryHDU(final_mask)  # we create a PrimaryHDU object to encapsulate the data
    wcs = WCS(header)
    new_hdu.header.update(wcs.to_header())
    new_hdu.writeto(i[:-5]+'_final.fits', overwrite=True) # write to a new file

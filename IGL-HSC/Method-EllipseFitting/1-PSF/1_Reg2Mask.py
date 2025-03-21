""""
=============================================================================
    Programa que genera una mascara bool a partir de las regiones 
                    DS9 'bright mask' de HSC
=============================================================================
"""


import matplotlib.pyplot as plt

from astropy.io import fits
from regions import read_ds9
from astropy.wcs import WCS
import numpy as np

# *****************************************************************************
# TO DO 
# *****************************************************************************
# To do: crear criterio que seleccione automaticamente esas regiones en All.reg
# ?? Actualizar paths a folder local de 'regions', 'writeto', 'data' y 'header': /Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/data/brightStars/S18A_starmask/new_S18Amask_g.reg
# ?? Hacer para cada banda individualmente

# =============================================================================
# Cargamos los datos 
# =============================================================================
results_dir = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'

data = fits.getdata('/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/data/images/400138-gri.fits')
shape = data.shape
header = fits.getheader('/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/data/images/400138-gri.fits')
wcs = WCS(header)

# =============================================================================
# Leemos todas las regiones
# =============================================================================
star_number = 4

regions = read_ds9(results_dir+'PSF/ds9Regions/star'+str(star_number)+'.reg')

# =============================================================================
# Hacemos un bucle para crear la máscara
# =============================================================================

mask = np.zeros_like(data)
for i in regions:    
   mask = mask +i.to_pixel(wcs).to_mask('center').to_image(shape)

mask = np.int_(np.bool_(mask))

plt.imshow(mask,origin='lower')
plt.show()
mask_hdu = fits.PrimaryHDU(mask)  # we create a PrimaryHDU object to encapsulate the data
wcs = WCS(header)
mask_hdu.header.update(wcs.to_header())
mask_hdu.writeto(results_dir+'PSF/masks/mask_s'+str(star_number)+'.fits', overwrite=True) # write to a new file

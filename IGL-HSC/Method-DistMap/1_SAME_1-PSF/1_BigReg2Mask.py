""""
=============================================================================
    Programa que genera una mascara bool a partir de regiones DS9
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
# Este codigo coge un archivo 1.reg donde solo estan las regiones de mi fov
# To do: crear criterio que seleccione automaticamente esas regiones en All.reg
# Actualizar paths a folder local de 'regions', 'writeto', 'data' y 'header': /Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/data/brightStars/S18A_starmask/new_S18Amask_g.reg
# Hacer para cada banda individualmente
# guardar mascaraas como brightObjs-band.fits en /Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/400138/masks/intermediate/
# guardar regiones seleccionadas

groupID = '400138'

PSF_dir = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/PSF/'


# =============================================================================
# Cargamos los datos 
# =============================================================================
mask = fits.getdata(PSF_dir+'masks/mask_s2_final_intermediate.fits')
shape = mask.shape
header = fits.getheader(PSF_dir+'masks/mask_s2_final_intermediate.fits')
wcs = WCS(header)

# =============================================================================
# Leemos todas las regiones
# =============================================================================

reg = read_ds9(PSF_dir+'ds9Regions/bigreg2.reg')

# =============================================================================
# Hacemos un bucle para crear la máscara
# =============================================================================

for i in reg:    
   mask = mask +i.to_pixel(wcs).to_mask('center').to_image(shape)
   
mask[mask!=0]=np.nan   
   

new_hdu = fits.PrimaryHDU(mask)  # we create a PrimaryHDU object to encapsulate the data
wcs = WCS(header)
new_hdu.header.update(wcs.to_header())
new_hdu.writeto(PSF_dir+'masks/mask_s2_final.fits', overwrite=True) # write to a new file

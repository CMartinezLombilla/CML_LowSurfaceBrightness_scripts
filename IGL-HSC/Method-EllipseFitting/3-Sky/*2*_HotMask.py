#*****************************************************************************
#                       COLD MASK IMAGE PROGRAM   [Rix+2004]
#-----------------------------------------------------------------------------
#                             Cristina Martinez 08/20
#*****************************************************************************

# Versions
# %pip install numpy== 1.18.5
# %pip install matplotlib==3.2.2
# %pip install astropy==4.0
# %pip install photutils==0.7.2

#!pip show numpy
#!pip show matplotlib
#!pip show astropy
#!pip show photutils


# Basics
import numpy as np

# Matplotlib
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
plt.close('all')
plt.ion()

# Astropy
from photutils import detect_threshold
from astropy.convolution import Gaussian2DKernel
from photutils import detect_sources
from photutils import deblend_sources
from astropy.convolution import convolve

# CML
import sys
#sys.path.append('/Users/felipe/Dropbox/CML/')
sys.path.append('/Users/z3530031/Dropbox/CML/')

from CML import DeleteMaskRegion
from CML import Patchs



# =============================================================================
# Esta lineas es para saber el tiempo de ejecución
# =============================================================================

import time
start = time.time()


# =============================================================================
# Cargamos las imagenes y las mascaras 
# =============================================================================

# Group info
groupID = '400138'
band = 'gri'  

# Paths & data
masks_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/masks/intermediate/'
input_mask = fits.getdata(masks_path+'unmaskedCore.fits')

image = fits.getdata('/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/PSF/i-star3.fits')
header = fits.getheader('/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/PSF/i-star3.fits')


# Convolucionamos imagen gri con kernel gaussiano [to creat an unsharp image]
sigma_pix = 5
kernel= Gaussian2DKernel(sigma_pix) #X & Y size by default 8*sigma_pix
conv_im = convolve(image, kernel)

# Creamos mascara bool 
bool_mask = np.isnan(input_mask)

# restamos gri - conv_gri
unsharp = image - conv_im

# enmascaramos con la mascara cold
data = np.zeros_like(unsharp)
data[bool_mask] = np.nan
data[~bool_mask] = unsharp[~bool_mask]


# =============================================================================
# Sorce detection & deblending en secciones mas pequenas sobre la imagen
# =============================================================================
patchs = Patchs(data,2100)

# Construimos una matriz de 0 del tamaño de la imagen como estructura para la mascara hot
mask = np.zeros(np.shape(data))

for i in patchs:
	# Establish a Hot nsigma threshold (we detect everything above)
	threshold = detect_threshold(data[i], nsigma=1.5, sigclip_sigma=3., mask_value = np.nan)
	# Segmentation - detecting sources **********************************************
	segm = detect_sources(data[i], threshold, npixels=9, mask=bool_mask[i])
	# Deblending sources **********************************************
	#mask[i] = deblend_sources(data[i], segm, npixels=7, filter_kernel=kernel).data
	mask[i] = deblend_sources(data[i], segm, npixels=9).data


# =============================================================================
# Fin contador de tiempo
# =============================================================================
print('It took', time.time()-start, 'seconds.')



# =============================================================================
# Quitamos máscaras de sobra
# =============================================================================
print('Remove masks over the BCG')

patchs = Patchs(data,2200) # Cambiamos el tamaño de los parches


for i in patchs:
    # Editamos a mano la seccion y la reemplazamos en la original
    mask[i] = DeleteMaskRegion(mask[i],data[i], vmax=1.5)




# =============================================================================
# Convolucionamos la mascara con un kernel Gaussiano para agrandarla 
# =============================================================================

sigma_pix = 1
kernel= Gaussian2DKernel(sigma_pix, x_size=3, y_size=3) 

conv_mask = convolve(mask, kernel)

# =============================================================================
# Enmascaramos la imagen con la mascara resultante
# =============================================================================

data[conv_mask != 0]= np.nan
input_mask[conv_mask != 0] = np.nan



# =============================================================================
# Guardamos la imagen enmascarada 
# =============================================================================

masked_hdu = fits.PrimaryHDU(input_mask)  # we create a PrimaryHDU object to encapsulate the data:
wcs = WCS(header)
masked_hdu.header.update(wcs.to_header())
masked_hdu.header.append(('HISTORY', 'Hot + Cold Mask', 'C. Martinez-Lombilla'), end=True)
masked_hdu.writeto(masks_path+'HotColdMasked.fits', overwrite=True) # write to a new file



# =============================================================================
# Guardamos la imagen unshaped enmascarada 
# =============================================================================


masked_hdu = fits.PrimaryHDU(data)  # we create a PrimaryHDU object to encapsulate the data:
wcs = WCS(header)
masked_hdu.header.update(wcs.to_header())
masked_hdu.header.append(('HISTORY', 'Hot + Cold Mask Ushaped', 'C. Martinez-Lombilla'), end=True)
masked_hdu.writeto(masks_path+'UnshapedMasked.fits', overwrite=True) # write to a new file















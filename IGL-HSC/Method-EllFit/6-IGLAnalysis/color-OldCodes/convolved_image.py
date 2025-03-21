
import numpy as np
import os 
from astropy.io import fits
import matplotlib.pyplot as plt
plt.ion()
import math
import matplotlib.pylab as plb


### Cambiar de carpeta en a terminal si quiero hacerlo en otro FILTRO/GALAXIA !!!!!!!


galaxy='UGC01040'
param = np.loadtxt('datos.txt')
Xc=param[16]
Yc=param[15]
DMpc=param[85]

scale = 0.39612 # escala del pixel de SDSS
zp =24.+2.5*np.log10(53.907456)+5*np.log10(scale) # zero point de la imagen de Stripe82

#counts_r=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/deconvolved.fits")   # modelo_sinPSF+residuals
counts_r=fits.open("new_image_grideep.fits")  # IMAGEN ENMASCARADA
#counts_g=fits.open("rotated.fits")   #SIN ENMASCARAR
copy_counts_r=counts_r[0].data[:]

un_pix_en_kpc=math.tan((scale/(60.*60.))*(math.pi/180.))*DMpc*1.e3
un_kpc_en_pix=1/un_pix_en_kpc
un_kpc_en_arcsec=un_kpc_en_pix*scale  # hallo a cuantos pixeles (y arcsec) equivale una escala fisica de resolucion (1 kpc)


##********************************
## RESOLUCION EQUIVALENTE A 1kpc**
##********************************

fwhm=un_kpc_en_pix
sigma_gauss=fwhm/2.35482

from astropy.convolution import Gaussian2DKernel, convolve
gauss = Gaussian2DKernel(stddev=sigma_gauss)
output_convolution_r = convolve(copy_counts_r, gauss, boundary='extent')  #NaN values are replaced with interpolated values using the kernel as an interpolation function.


n=np.arange(len(counts_r[0].data[:,0])*1.*len(counts_r[0].data[0,:])*1.).reshape(len(counts_r[0].data[:,0]),len(counts_r[0].data[0,:]))
#hdu = fits.PrimaryHDU(n)
#hdu.writeto('convolved_1kpc.fits')
gauss_convolution_r=fits.open("convolved_1kpc.fits", mode='update')
gauss_convolution_r[0].header= counts_r[0].header
gauss_convolution_r[0].header['HISTORY'] = 'Convolved gaussian kernel 1kpc'

gauss_convolution_r[0].data[:]=output_convolution_r[:]
gauss_convolution_r.flush()





##********************************
## RESOLUCION EQUIVALENTE A 0.5kpc**
##********************************

half_kpc_en_pix=un_kpc_en_pix/2.

fwhm=half_kpc_en_pix
sigma_gauss=fwhm/2.35482

from astropy.convolution import Gaussian2DKernel, convolve
gauss = Gaussian2DKernel(stddev=sigma_gauss)
output_convolution_r = convolve(copy_counts_r, gauss, boundary='extent')  #NaN values are replaced with interpolated values using the kernel as an interpolation function.

n=np.arange(len(counts_r[0].data[:,0])*1.*len(counts_r[0].data[0,:])*1.).reshape(len(counts_r[0].data[:,0]),len(counts_r[0].data[0,:]))
#hdu = fits.PrimaryHDU(n)
#hdu.writeto('convolved_0.5kpc.fits')
counts_r=fits.open("new_image_grideep.fits")  
gauss_convolution_r=fits.open("convolved_0.5kpc.fits", mode='update')
gauss_convolution_r[0].header= counts_r[0].header
gauss_convolution_r[0].header['HISTORY'] = 'Convolved gaussian kernel 0.5kpc'

gauss_convolution_r[0].data[:]=output_convolution_r[:]
gauss_convolution_r.flush()



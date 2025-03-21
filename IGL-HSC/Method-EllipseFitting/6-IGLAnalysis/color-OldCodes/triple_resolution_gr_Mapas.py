
#####
#			************
#
#   CAMBIAR RUTAS DE ARCHIVOS AL CAMBIAR DE GALAXIA / FILTROS !!!!!!!!!!!
#
#			************
#####


import numpy as np
import os 
from astropy.io import fits
import matplotlib.pyplot as plt
plt.ion()
import math
import matplotlib.pylab as plb

galaxy='UGC01040'
param = np.loadtxt('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/datos.txt')
Xc=param[16]
Yc=param[15]
DMpc=param[85]

scale = 0.39612 # escala del pixel de SDSS
zp =24.+2.5*np.log10(53.907456)+5*np.log10(scale) # zero point de la imagen de Stripe82



#DATOS CONVOLUCIONADOS (cuentas)
################################
################################

gauss_convolution_1r=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/convolved_1kpc.fits")  # IMAGEN ENMASCARADA
gauss_convolution_1g=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/convolved_1kpc.fits")  # IMAGEN ENMASCARADA

gauss_convolution_05r=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/convolved_0.5kpc.fits")  # IMAGEN ENMASCARADA
gauss_convolution_05g=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/convolved_0.5kpc.fits")  # IMAGEN ENMASCARADA


#DATOS EN RESOLUCION ORIGINAL  (mag/arcsec2 y cuentas)
######################################################
######################################################

#counts_r=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/deconvolved.fits")   # modelo_sinPSF+residuals
#counts_g=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/deconvolved.fits")   # modelo_sinPSF+residuals

counts_r=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/new_image_grideep.fits")  # IMAGEN ENMASCARADA
counts_g=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/new_image_grideep.fits")  # IMAGEN ENMASCARADA

#counts_r=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/rotated.fits")   #SIN ENMASCARAR
#counts_g=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/rotated.fits")   #SIN ENMASCARAR


mu_r=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/mapa_MagArcsec2.fits")
mu_g=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/mapa_MagArcsec2.fits") #Voy haciendo todos en g y en r porque lo necesito para los mapas de color/masa/...



# IMAGEN EN TRIPLE RESOLUCION (en mag/arcsec2 y en cuentas)
##########################################################
##########################################################

#Primero genero la imagen entre 25 - inf magnitudes/equiv.en cuentas:

n=np.arange(len(mu_r[0].data[:,0])*1.*len(mu_r[0].data[0,:])*1.).reshape(len(mu_r[0].data[:,0]),len(mu_r[0].data[0,:]))

hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/counts_triple_resolution_gr.fits')
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/mu_triple_resolution_gr.fits')
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/counts_triple_resolution_gr.fits')
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/mu_triple_resolution_gr.fits')

counts_triple_resolution_r=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/counts_triple_resolution_gr.fits', mode='update')
mu_triple_resolution_r=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/mu_triple_resolution_gr.fits', mode='update')
counts_triple_resolution_g=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/counts_triple_resolution_gr.fits', mode='update')
mu_triple_resolution_g=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/mu_triple_resolution_gr.fits', mode='update')


for x in range(len(mu_r[0].data[:,0])):  
  for y in range(len(mu_r[0].data[0,:])):
	if -2.5*(np.log10(gauss_convolution_1r[0].data[x,y]))+zp>=25.5:
		mu_triple_resolution_r[0].data[x,y]=-2.5*(np.log10(gauss_convolution_1r[0].data[x,y]))+zp
		counts_triple_resolution_r[0].data[x,y]=gauss_convolution_1r[0].data[x,y]
		mu_triple_resolution_g[0].data[x,y]=-2.5*(np.log10(gauss_convolution_1g[0].data[x,y]))+zp
		counts_triple_resolution_g[0].data[x,y]=gauss_convolution_1g[0].data[x,y]
	elif 25.5<-2.5*(np.log10(gauss_convolution_05r[0].data[x,y]))+zp>=24.5:
		mu_triple_resolution_r[0].data[x,y]=-2.5*(np.log10(gauss_convolution_05r[0].data[x,y]))+zp
		counts_triple_resolution_r[0].data[x,y]=gauss_convolution_05r[0].data[x,y]
		mu_triple_resolution_g[0].data[x,y]=-2.5*(np.log10(gauss_convolution_05g[0].data[x,y]))+zp
		counts_triple_resolution_g[0].data[x,y]=gauss_convolution_05g[0].data[x,y]	
	else:
		mu_triple_resolution_r[0].data[x,y]=mu_r[0].data[x,y]
		counts_triple_resolution_r[0].data[x,y]=counts_r[0].data[x,y]
		mu_triple_resolution_g[0].data[x,y]=mu_g[0].data[x,y]
		counts_triple_resolution_g[0].data[x,y]=counts_g[0].data[x,y]


mu_triple_resolution_r[0].header=mu_r[0].header
mu_triple_resolution_r[0].header['HISTORY'] = 'triple resolution [mag/arcsec2] (original + 0.5kpc/pix + 1kpc/pix)'
mu_triple_resolution_r.flush()
counts_triple_resolution_r[0].header=counts_r[0].header
counts_triple_resolution_r[0].header['HISTORY'] = 'triple resolution [counts] (original + 0.5kpc/pix + 1kpc/pix)'
counts_triple_resolution_r.flush()
mu_triple_resolution_g[0].header=mu_g[0].header
mu_triple_resolution_g[0].header['HISTORY'] = 'triple resolution [mag/arcsec2]  (original + 0.5kpc/pix + 1kpc/pix)'
mu_triple_resolution_g.flush()
counts_triple_resolution_g[0].header=counts_g[0].header
counts_triple_resolution_g[0].header['HISTORY'] = 'triple resolution [counts] (original + 0.5kpc/pix + 1kpc/pix)'
counts_triple_resolution_g.flush()




#plt.imshow(copy_mu)
#plt.figure(1)
#plt.pcolor(copy_mu, vmin=np.nanmax(copy_mu), vmax=np.nanmin(copy_mu))
#plt.title('triple resolution Mag/arcsec$^{2}$ map of %s' %galaxy)
#plt.colorbar()


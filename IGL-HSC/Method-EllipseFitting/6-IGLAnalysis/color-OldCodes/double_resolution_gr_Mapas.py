
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

gauss_convolution_r=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/convolved_1kpc.fits")  # IMAGEN ENMASCARADA
gauss_convolution_g=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/convolved_1kpc.fits")  # IMAGEN ENMASCARADA


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



# IMAGEN EN DOBLE RESOLUCION (en mag/arcsec2 y en cuentas)
##########################################################
##########################################################

#Primero genero la imagen entre 25 - inf magnitudes/equiv.en cuentas:

copy_mu_r=np.arange(len(mu_r[0].data[:,0])*1.*len(mu_r[0].data[0,:])*1.).reshape(len(mu_r[0].data[:,0]),len(mu_r[0].data[0,:]))
copy_counts_r=np.arange(len(counts_r[0].data[:,0])*1.*len(counts_r[0].data[0,:])*1.).reshape(len(counts_r[0].data[:,0]),len(counts_r[0].data[0,:]))

copy_mu_g=np.arange(len(mu_g[0].data[:,0])*1.*len(mu_g[0].data[0,:])*1.).reshape(len(mu_g[0].data[:,0]),len(mu_g[0].data[0,:]))
copy_counts_g=np.arange(len(counts_g[0].data[:,0])*1.*len(counts_g[0].data[0,:])*1.).reshape(len(counts_g[0].data[:,0]),len(counts_g[0].data[0,:]))

for x in range(len(mu_r[0].data[:,0])):  
  for y in range(len(mu_r[0].data[0,:])):
	if -2.5*(np.log10(gauss_convolution_r[0].data[x,y]))+zp>=25.0:
		copy_mu_r[x,y]=-2.5*(np.log10(gauss_convolution_r[0].data[x,y]))+zp
		copy_counts_r[x,y]=gauss_convolution_r[0].data[x,y]
		copy_mu_g[x,y]=-2.5*(np.log10(gauss_convolution_g[0].data[x,y]))+zp
		copy_counts_g[x,y]=gauss_convolution_g[0].data[x,y]
	else:
		copy_mu_r[x,y]=0.
		copy_counts_r[x,y]=0.
		copy_mu_g[x,y]=0.
		copy_counts_g[x,y]=0.

#plt.imshow(copy_mu)
#plt.figure(1)
#plt.pcolor(copy_mu, vmin=np.nanmax(copy_mu), vmax=np.nanmin(copy_mu))
#plt.title('Double resolution Mag/arcsec$^{2}$ map of %s' %galaxy)
#plt.colorbar()

n=np.arange(len(mu_r[0].data[:,0])*1.*len(mu_r[0].data[0,:])*1.).reshape(len(mu_r[0].data[:,0]),len(mu_r[0].data[0,:]))
hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/mu_major25.0_gr.fits')
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/counts_major25.0_gr.fits')
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/mu_major25.0_gr.fits')
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/counts_major25.0_gr.fits')
mu_major25_r=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/mu_major25.0_gr.fits', mode='update')
counts_major25_r=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/counts_major25.0_gr.fits', mode='update')
mu_major25_g=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/mu_major25.0_gr.fits', mode='update')
counts_major25_g=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/counts_major25.0_gr.fits', mode='update')

mu_major25_r[0].data[:]=copy_mu_r[:]
gauss_convolution_r=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/convolved_1kpc.fits")
mu_major25_r[0].header=gauss_convolution_r[0].header
mu_major25_r[0].header['HISTORY'] = 'surface brightness great than 25'
mu_major25_r.flush()
counts_major25_r[0].data[:]=copy_counts_r[:]
gauss_convolution_r=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/convolved_1kpc.fits")
counts_major25_r[0].header= gauss_convolution_r[0].header
counts_major25_r[0].header['HISTORY'] = 'surface brightness great than 25'
counts_major25_r.flush()
mu_major25_g[0].data[:]=copy_mu_g[:]
gauss_convolution_g=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/convolved_1kpc.fits")
mu_major25_g[0].header= gauss_convolution_g[0].header
mu_major25_g[0].header['HISTORY'] = 'surface brightness great than 25'
mu_major25_g.flush()
counts_major25_g[0].data[:]=copy_counts_g[:]
gauss_convolution_g=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/convolved_1kpc.fits")
counts_major25_g[0].header= gauss_convolution_g[0].header
counts_major25_g[0].header['HISTORY'] = 'surface brightness great than 25'
counts_major25_g.flush()


# Ahora creo la imagen completa a doble resolucion

mu_major25_r=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/mu_major25.0_gr.fits')
mu_major25_g=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/mu_major25.0_gr.fits')

counts_major25_r=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/counts_major25.0_gr.fits')
counts_major25_g=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/counts_major25.0_gr.fits')

gauss_convolution_r=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/convolved_1kpc.fits")
gauss_convolution_g=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/convolved_1kpc.fits")

n=np.arange(len(mu_major25_r[0].data[:,0])*1.*len(mu_major25_r[0].data[0,:])*1.).reshape(len(mu_major25_r[0].data[:,0]),len(mu_major25_r[0].data[0,:]))
hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/counts_double_resolution_gr.fits')
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/mu_double_resolution_gr.fits')
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/counts_double_resolution_gr.fits')
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/mu_double_resolution_gr.fits')

counts_double_resolution_r=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/counts_double_resolution_gr.fits', mode='update')
mu_double_resolution_r=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/mu_double_resolution_gr.fits', mode='update')
counts_double_resolution_g=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/counts_double_resolution_gr.fits', mode='update')
mu_double_resolution_g=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/mu_double_resolution_gr.fits', mode='update')

for x in range(len(mu_major25_r[0].data[:,0])):  
  for y in range(len(mu_major25_r[0].data[0,:])):
	if mu_major25_r[0].data[x,y]==0.:
		mu_double_resolution_r[0].data[x,y]=mu_r[0].data[x,y]
		mu_double_resolution_g[0].data[x,y]=mu_g[0].data[x,y]
		counts_double_resolution_r[0].data[x,y]=counts_r[0].data[x,y]
		counts_double_resolution_g[0].data[x,y]=counts_g[0].data[x,y]		

	else:
		mu_double_resolution_r[0].data[x,y]=mu_major25_r[0].data[x,y]
		mu_double_resolution_g[0].data[x,y]=mu_major25_g[0].data[x,y]
		counts_double_resolution_r[0].data[x,y]=counts_major25_r[0].data[x,y]
		counts_double_resolution_g[0].data[x,y]=counts_major25_g[0].data[x,y]

mu_double_resolution_r[0].header=gauss_convolution_r[0].header
mu_double_resolution_r[0].header['HISTORY'] = 'double resolution (original + 1kpc/pix)'
mu_double_resolution_r.flush()
counts_double_resolution_r[0].header=gauss_convolution_r[0].header
counts_double_resolution_r[0].header['HISTORY'] = 'double resolution (original + 1kpc/pix)'
counts_double_resolution_r.flush()
mu_double_resolution_g[0].header=gauss_convolution_r[0].header
mu_double_resolution_g[0].header['HISTORY'] = 'double resolution (original + 1kpc/pix)'
mu_double_resolution_g.flush()
counts_double_resolution_g[0].header=gauss_convolution_r[0].header
counts_double_resolution_g[0].header['HISTORY'] = 'double resolution (original + 1kpc/pix)'
counts_double_resolution_g.flush()



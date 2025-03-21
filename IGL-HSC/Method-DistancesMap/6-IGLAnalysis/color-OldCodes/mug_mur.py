import numpy as np
import os 
from astropy.io import fits
import matplotlib.pyplot as plt
plt.ion()
import math
import matplotlib.pylab as plb

galaxy='UGC01040'
param = np.loadtxt('datos.txt')
DMpc=param[85]
bulge_lim=param[86]
semi_mayor=param[87]
semi_menor=param[88]
log_d25=param[89]
Xc=param[16]
Yc=param[15]

### RESTAMOS LOS MAPAS DE MAGNITUD PARA COMPROBAR QUE EL COLOR SALE IGUAL QUE CON LA OTRA FORMULITA DE LOS LOGARITMOS
scale = 0.39612 # escala del pixel de SDSS
zp =24.+2.5*np.log10(53.907456)+5*np.log10(scale) # zero point de la imagen de Stripe82


mu_double_r=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/mu_double_resolution_r.fits', mode='update')
mu_double_g=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/mu_double_resolution_g.fits', mode='update')

color_g_r_double=np.arange(len(mu_double_r[0].data[:,0])*1.*len(mu_double_r[0].data[0,:])*1.).reshape(len(mu_double_r[0].data[:,0]),len(mu_double_r[0].data[0,:]))

for x in range(len(mu_double_r[0].data[:,0])):  
  for y in range(len(mu_double_r[0].data[0,:])):
	color_g_r_double[x,y]=mu_double_g[0].data[x,y]-mu_double_r[0].data[x,y]

n=np.arange(len(mu_double_r[0].data[:,0])*1.*len(mu_double_r[0].data[0,:])*1.).reshape(len(mu_double_r[0].data[:,0]),len(mu_double_r[0].data[0,:]))
hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/prueba_color.fits')
mapa_color_g_r_double=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/prueba_color.fits', mode='update')
mapa_color_g_r_double[0].data[:]=color_g_r_double[:]
mapa_color_g_r_double.flush()



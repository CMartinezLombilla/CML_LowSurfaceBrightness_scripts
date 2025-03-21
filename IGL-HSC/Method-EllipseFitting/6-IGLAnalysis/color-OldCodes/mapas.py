#####
#			************
#
#   CAMBIAR RUTAS DE ARCHIVOS AL CAMBIAR DE GALAXIA / FILTROS !!!!!!!!!!!
#		PONERME EN LA CARPETA DEL FILTRO
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
param = np.loadtxt('datos.txt')
DMpc=param[85]
bulge_lim=param[86]
semi_mayor=param[87]
semi_menor=param[88]
log_d25=param[89]
Xc=param[16]
Yc=param[15]


scale = 0.39612 # escala del pixel de SDSS
zp =24.+2.5*np.log10(53.907456)+5*np.log10(scale) # zero point de la imagen de Stripe82



#GENERO EL MAPA DE MAGNITUDES [MAG/ARCSEC2]
#******************************************

#DATOS

#image=fits.open("deconvolved.fits")   # modelo_sinPSF+residuals !!!
image=fits.open("new_image_grideep.fits")  # IMAGEN ENMASCARADA
#image=fits.open("rotated.fits")   # SIN ENMASCARAR !!!
copyimage_mag=np.arange(len(image[0].data[:,0])*1.*len(image[0].data[0,:])*1.).reshape(len(image[0].data[:,0]),len(image[0].data[0,:]))
for x in range(len(image[0].data[:,0])):  
  for y in range(len(image[0].data[0,:])):
	copyimage_mag[x,y]=-2.5*(np.log10(image[0].data[x,y]))+zp

mapa_magimage=copyimage_mag

#plt.figure(1)
#plt.imshow(mapa_magimage)
#plt.pcolor(mapa_magimage, vmin=np.nanmax(mapa_magimage), vmax=np.nanmin(mapa_magimage))
#plt.title('Mag/arcsec$^{2}$ map of %s' %galaxy)
#plt.colorbar()

n=np.arange(len(image[0].data[:,0])*1.*len(image[0].data[0,:])*1.).reshape(len(image[0].data[:,0]),len(image[0].data[0,:]))
hdu = fits.PrimaryHDU(n)
hdu.writeto('mapa_MagArcsec2.fits')
mapa=fits.open("mapa_MagArcsec2.fits", mode='update')
mapa[0].data[:]=mapa_magimage[:]
mapa[0].header= image[0].header
mapa[0].header['HISTORY'] = 'Surface brightness map [mag/arcsec2]'
mapa.flush()

################


modelPSF=fits.open("modelimagepsf.fits")
model=fits.open("modelimage.fits")


#MODELO con PSF

copyPSF_mag=np.arange(len(modelPSF[0].data[:,0])*1.*len(modelPSF[0].data[0,:])*1.).reshape(len(modelPSF[0].data[:,0]),len(modelPSF[0].data[0,:]))
for x in range(len(modelPSF[0].data[:,0])):  
  for y in range(len(modelPSF[0].data[0,:])):
	copyPSF_mag[x,y]=-2.5*np.log10(modelPSF[0].data[x,y])+zp

mapa_magPSF=copyPSF_mag


#MODELO

copy_mag=np.arange(len(model[0].data[:,0])*1.*len(model[0].data[0,:])*1.).reshape(len(model[0].data[:,0]),len(model[0].data[0,:]))
for x in range(len(model[0].data[:,0])):  
  for y in range(len(model[0].data[0,:])):
	copy_mag[x,y]=-2.5*np.log10(model[0].data[x,y])+zp
	
mapa_mag=copy_mag



plt.figure(2)
plt.imshow(mapa_magPSF)
plt.pcolor(mapa_magPSF, vmin=np.max(mapa_magPSF), vmax=np.min(mapa_magPSF))
plt.title('Mag/arcsec$^{2}$ map of %s (PSF convolved model)' %galaxy)
plt.colorbar()
plt.figure(3)
plt.imshow(mapa_mag)
plt.pcolor(mapa_mag, vmin=np.max(mapa_magPSF), vmax=np.min(mapa_magPSF))
plt.title('Mag/arcsec$^{2}$ map of %s (Deconvolved model)' %galaxy)
plt.colorbar()



#GENERO EL MAPA DE DISTANCIAS [KPC DESDE EL PLANO]
#*************************************************

Dkpc=DMpc*1e3


#MODELO con PSF

copyPSF_kpc=np.arange(len(modelPSF[0].data[:,0])*1.*len(modelPSF[0].data[0,:])*1.).reshape(len(modelPSF[0].data[:,0]),len(modelPSF[0].data[0,:]))
for x in range(len(modelPSF[0].data[:,0])):  
  for y in range(len(modelPSF[0].data[0,:])):
	to_arcsec=abs(x-Xc)*scale
	to_rad=(to_arcsec/(60.*60.))*(math.pi/180.)
	x_kpc=math.tan(to_rad)*Dkpc
	copyPSF_kpc[x,y]=x_kpc

mapa_kpcPSF=copyPSF_kpc



#MODELO 

copy_kpc=np.arange(len(model[0].data[:,0])*1.*len(model[0].data[0,:])*1.).reshape(len(model[0].data[:,0]),len(model[0].data[0,:]))
for x in range(len(model[0].data[:,0])):  
  for y in range(len(model[0].data[0,:])):
	to_arcsec=abs(x-Xc)*scale
	to_rad=(to_arcsec/(60.*60.))*(math.pi/180.)
	x_kpc=math.tan(to_rad)*Dkpc
	copy_kpc[x,y]=x_kpc

mapa_kpc=copy_kpc  #lo logico es que ambos mapas de distancia sean iguales, porque ambas imagenes tienen el mismo tamano


plt.figure(4)
plt.imshow(mapa_kpcPSF)
plt.title('Distance [kpc] from the midplane map of %s \n (PSF convolved model)' %galaxy)
plt.colorbar()
plt.figure(5)
plt.imshow(mapa_kpc)
plt.title('Distance [kpc] from the midplane map of %s \n (Deconvolved model)' %galaxy)
plt.colorbar()



#***********
#CALCULOS
#***********


#FLUJO TOTAL DE LA GALAXIA
#******************************************

mag_lim=30.

l=[]
for x in range(len(mapa_magPSF[:,0])):  
  for y in range(len(mapa_magPSF[0,:])):
	if mapa_magPSF[x,y]<mag_lim:
		l.append(modelPSF[0].data[x,y])
total_fluxPSF=sum(l)

s=[]
for x in range(len(mapa_mag[:,0])):  
  for y in range(len(mapa_mag[0,:])):
	if mapa_mag[x,y]<mag_lim:
		s.append(model[0].data[x,y])
total_flux=sum(s)


#FLUJO DE LA GALAXIA > 1kpc
#******************************************

list_lim=[1., bulge_lim, 2., 3., 4.]

for kpc_lim in list_lim:
	m=[]
	for x in range(len(mapa_magPSF[:,0])):    
		for y in range(len(mapa_magPSF[0,:])):
			if mapa_kpcPSF[x,y]>kpc_lim and mapa_magPSF[x,y]<mag_lim:
				m.append(modelPSF[0].data[x,y])
	kpc_fluxPSF=sum(m)

	t=[]
	for x in range(len(mapa_mag[:,0])):    
		for y in range(len(mapa_mag[0,:])):
			if mapa_kpc[x,y]>kpc_lim and mapa_mag[x,y]<mag_lim:
				t.append(model[0].data[x,y])
	kpc_flux=sum(t)

	#COCIENTES DE FLUJOS
	#******************************************

	ratioPSF=kpc_fluxPSF/total_fluxPSF*100.
	ratio=kpc_flux/total_flux*100.
	
	print 'Distancia limite: '+str(kpc_lim)+' kpc'
	print 'En el modelo convolucionado con la PSF, el '+str(ratioPSF)+'% de la luz se concentra en las partes externas'
	print 'En el modelo deconcolucionado, el '+str(ratio)+'% de la luz se concentra en las partes externas'



#GENERO EL MAPA DE MASAS [kg]  (g-r)
#******************************************

# 1.  obtain M/L ratio from the (g-r) colors + stellar density

#counts_r=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/deconvolved.fits")   # modelo_sinPSF+residuals
#counts_g=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/deconvolved.fits")   # modelo_sinPSF+residuals

counts_r=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/new_image_grideep.fits")  # IMAGEN ENMASCARADA
counts_g=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/new_image_grideep.fits")  # IMAGEN ENMASCARADA

a_lambda=-0.306
b_lambda=1.097
m_sun_r=4.76


log_M_to_L=np.arange(len(counts_r[0].data[:,0])*1.*len(counts_r[0].data[0,:])*1.).reshape(len(counts_r[0].data[:,0]),len(counts_r[0].data[0,:]))
M_to_L=np.arange(len(counts_r[0].data[:,0])*1.*len(counts_r[0].data[0,:])*1.).reshape(len(counts_r[0].data[:,0]),len(counts_r[0].data[0,:]))

color_g_r=np.arange(len(counts_r[0].data[:,0])*1.*len(counts_r[0].data[0,:])*1.).reshape(len(counts_r[0].data[:,0]),len(counts_r[0].data[0,:]))

log_stelar_mass_desity=np.arange(len(counts_r[0].data[:,0])*1.*len(counts_r[0].data[0,:])*1.).reshape(len(counts_r[0].data[:,0]),len(counts_r[0].data[0,:]))


for x in range(len(counts_r[0].data[:,0])):  
  for y in range(len(counts_r[0].data[0,:])):
	color_g_r[x,y]=-2.5*np.log10(counts_g[0].data[x,y]/counts_r[0].data[x,y]) # mismo zp para todos los filtros en Stripe82
	log_M_to_L[x,y]=(a_lambda+b_lambda*color_g_r[x,y])-0.15
	log_stelar_mass_desity[x,y]=log_M_to_L[x,y]-0.4*(-2.5*np.log10(counts_r[0].data[x,y])-m_sun_r)+8.629
"""		
plt.figure(6)
plt.imshow(color_g_r)
plt.pcolor(color_g_r,  vmin=np.nanmean(color_g_r)-3*np.nanstd(color_g_r),vmax=np.nanmax(color_g_r))
plt.title('%s: color $(g-r)_{ABmag}$' %galaxy)
plt.colorbar()
	
plt.figure(7)
plt.imshow(log_M_to_L)
plt.pcolor(log_M_to_L, vmin=np.nanmean(log_M_to_L)-3*np.nanstd(log_M_to_L),vmax=np.nanmax(log_M_to_L))
plt.title('%s: $log_{10}(M/L$) in r-band' %galaxy)
plt.colorbar()

plt.figure(8)
plt.imshow(log_stelar_mass_desity)
plt.pcolor(log_stelar_mass_desity, vmin=0.,vmax=np.nanmax(log_stelar_mass_desity[Xc-10.:Xc+10.,Yc-10.:Yc+10.]))
plt.title('%s: $log_{10}\Sigma$ [$M_{\odot} pc^{-2}$] in r-band' %galaxy)
plt.colorbar()
"""

n=np.arange(len(counts_r[0].data[:,0])*1.*len(counts_r[0].data[0,:])*1.).reshape(len(counts_r[0].data[:,0]),len(counts_r[0].data[0,:]))
hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_color_g_r.fits')
mapa_color_g_r=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_color_g_r.fits', mode='update')
mapa_color_g_r[0].data[:]=color_g_r[:]
mapa_color_g_r[0].header=counts_r[0].header
mapa_color_g_r[0].header['HISTORY'] = '(g-r) color map'
mapa_color_g_r.flush()


hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_log_M_to_L_r.fits')
mapa_log_M_to_L=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_log_M_to_L_r.fits', mode='update')
mapa_log_M_to_L[0].data[:]=log_M_to_L[:]
mapa_log_M_to_L[0].header=counts_r[0].header
mapa_log_M_to_L[0].header['HISTORY'] = 'r-band M/L map'
mapa_log_M_to_L.flush()


hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_StellarMassDensity_r.fits')
mapa_stelar_mass_desity=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_StellarMassDensity_r.fits', mode='update')
mapa_stelar_mass_desity[0].data[:]=log_stelar_mass_desity[:]
mapa_stelar_mass_desity[0].header=counts_r[0].header
mapa_stelar_mass_desity[0].header['HISTORY'] = 'r-band stellar mass density map'
mapa_stelar_mass_desity.flush()


#GENERO EL MAPA DE MASAS [kg]  (r-i)
#******************************************

# 1.  obtain M/L ratio from the (r-i) colors + stellar density

#counts_i=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/i/deconvolved.fits")   # modelo_sinPSF+residuals
#counts_r=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/deconvolved.fits")   # modelo_sinPSF+residuals

counts_i=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/i/new_image_grideep.fits")  # IMAGEN ENMASCARADA
counts_r=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/new_image_grideep.fits")  # IMAGEN ENMASCARADA

mu_i=fits.open("/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/i/mapa_MagArcsec2.fits")

a_lambda=-0.122
b_lambda=1.431
m_sun_i=4.58

log_M_to_L=np.arange(len(counts_i[0].data[:,0])*1.*len(counts_i[0].data[0,:])*1.).reshape(len(counts_i[0].data[:,0]),len(counts_i[0].data[0,:]))
M_to_L=np.arange(len(counts_i[0].data[:,0])*1.*len(counts_i[0].data[0,:])*1.).reshape(len(counts_i[0].data[:,0]),len(counts_i[0].data[0,:]))

color_r_i=np.arange(len(counts_i[0].data[:,0])*1.*len(counts_i[0].data[0,:])*1.).reshape(len(counts_i[0].data[:,0]),len(counts_i[0].data[0,:]))

log_stelar_mass_desity=np.arange(len(counts_i[0].data[:,0])*1.*len(counts_i[0].data[0,:])*1.).reshape(len(counts_i[0].data[:,0]),len(counts_i[0].data[0,:]))


for x in range(len(counts_i[0].data[:,0])):  
  for y in range(len(counts_i[0].data[0,:])):
	color_r_i[x,y]=-2.5*np.log10(counts_r[0].data[x,y]/counts_i[0].data[x,y]) # mismo zp para todos los filtros en Stripe82
	log_M_to_L[x,y]=(a_lambda+b_lambda*color_r_i[x,y])-0.15
	log_stelar_mass_desity[x,y]=log_M_to_L[x,y]-0.4*(-2.5*np.log10(counts_i[0].data[x,y])-m_sun_i)+8.629
"""		
plt.figure(9)
plt.imshow(color_r_i)
plt.pcolor(color_r_i,  vmin=np.nanmean(color_r_i)-3*np.nanstd(color_r_i),vmax=np.nanmax(color_r_i))
plt.title('%s: color $(g-r)_{ABmag}$' %galaxy)
plt.colorbar()
	
plt.figure(10)
plt.imshow(log_M_to_L)
plt.pcolor(log_M_to_L, vmin=np.nanmean(log_M_to_L)-3*np.nanstd(log_M_to_L),vmax=np.nanmax(log_M_to_L))
plt.title('%s: $log_{10}(M/L$) in i-band' %galaxy)
plt.colorbar()

plt.figure(11)
plt.imshow(log_stelar_mass_desity)
plt.pcolor(log_stelar_mass_desity, vmin=0.,vmax=np.nanmax(log_stelar_mass_desity[Xc-10.:Xc+10.,Yc-10.:Yc+10.]))
plt.title('%s: $log_{10}\Sigma$ [$M_{\odot} pc^{-2}$] in i-band' %galaxy)
plt.colorbar()
"""

n=np.arange(len(mu_i[0].data[:,0])*1.*len(mu_i[0].data[0,:])*1.).reshape(len(mu_i[0].data[:,0]),len(mu_i[0].data[0,:]))
hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_color_r_i.fits')
mapa_color_r_i=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_color_r_i.fits', mode='update')
mapa_color_r_i[0].data[:]=color_r_i[:]
mapa_color_r_i[0].header=counts_i[0].header
mapa_color_r_i[0].header['HISTORY'] = '(r-i) color map'
mapa_color_r_i.flush()


hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_log_M_to_L_i.fits')
mapa_log_M_to_L=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_log_M_to_L_i.fits', mode='update')
mapa_log_M_to_L[0].data[:]=log_M_to_L[:]
mapa_log_M_to_L[0].header=counts_i[0].header
mapa_log_M_to_L[0].header['HISTORY'] = 'i-band M/L map'
mapa_log_M_to_L.flush()


hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_StellarMassDensity_i.fits')
mapa_stelar_mass_desity=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_StellarMassDensity_i.fits', mode='update')
mapa_stelar_mass_desity[0].data[:]=log_stelar_mass_desity[:]
mapa_stelar_mass_desity[0].header=counts_i[0].header
mapa_stelar_mass_desity[0].header['HISTORY'] = 'i-band stellar mass density map'
mapa_stelar_mass_desity.flush()


##############################################################################################################
##############################################################################################################


##DOUBLE RESOLUTION  (mismo proceso pero a resolucion doble)
######################

# PARA EL COLOR (g-r)
counts_double_r=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/counts_double_resolution_gr.fits', mode='update')
mu_double_r=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/mu_double_resolution_gr.fits', mode='update')
counts_double_g=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/counts_double_resolution_gr.fits', mode='update')

a_lambda=-0.306
b_lambda=1.097
m_sun_r=4.76


log_M_to_L_double=np.arange(len(counts_double_r[0].data[:,0])*1.*len(counts_double_r[0].data[0,:])*1.).reshape(len(counts_double_r[0].data[:,0]),len(counts_double_r[0].data[0,:]))
M_to_L_double=np.arange(len(counts_double_r[0].data[:,0])*1.*len(counts_double_r[0].data[0,:])*1.).reshape(len(counts_double_r[0].data[:,0]),len(counts_double_r[0].data[0,:]))

color_g_r_double=np.arange(len(counts_double_r[0].data[:,0])*1.*len(counts_double_r[0].data[0,:])*1.).reshape(len(counts_double_r[0].data[:,0]),len(counts_double_r[0].data[0,:]))

log_stelar_mass_desity_double=np.arange(len(counts_double_r[0].data[:,0])*1.*len(counts_double_r[0].data[0,:])*1.).reshape(len(counts_double_r[0].data[:,0]),len(counts_double_r[0].data[0,:]))


for x in range(len(counts_double_r[0].data[:,0])):  
  for y in range(len(counts_double_r[0].data[0,:])):
	color_g_r_double[x,y]=-2.5*np.log10(counts_double_g[0].data[x,y]/counts_double_r[0].data[x,y]) # mismo zp para todos los filtros en Stripe82
	log_M_to_L_double[x,y]=(a_lambda+b_lambda*color_g_r_double[x,y])-0.15
	log_stelar_mass_desity_double[x,y]=log_M_to_L_double[x,y]-0.4*(-2.5*np.log10(counts_double_r[0].data[x,y])-m_sun_r)+8.629
		
n=np.arange(len(mu_double_r[0].data[:,0])*1.*len(mu_double_r[0].data[0,:])*1.).reshape(len(mu_double_r[0].data[:,0]),len(mu_double_r[0].data[0,:]))
hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_color_g_r_double.fits')
mapa_color_g_r_double=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_color_g_r_double.fits', mode='update')
mapa_color_g_r_double[0].data[:]=color_g_r_double[:]
mapa_color_g_r_double[0].header=counts_double_r[0].header
mapa_color_g_r_double[0].header['HISTORY'] = '(g-r) color map (double resolution)'
mapa_color_g_r_double.flush()


hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_log_M_to_L_r_double.fits')
mapa_log_M_to_L_double=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_log_M_to_L_r_double.fits', mode='update')
mapa_log_M_to_L_double[0].data[:]=log_M_to_L_double[:]
mapa_log_M_to_L_double[0].header=counts_double_r[0].header
mapa_log_M_to_L_double[0].header['HISTORY'] = 'r-band M/L map (double resolution)'
mapa_log_M_to_L_double.flush()


hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_StellarMassDensity_r_double.fits')
mapa_stelar_mass_desity_double=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_StellarMassDensity_r_double.fits', mode='update')
mapa_stelar_mass_desity_double[0].data[:]=log_stelar_mass_desity_double[:]
mapa_stelar_mass_desity_double[0].header=counts_double_r[0].header
mapa_stelar_mass_desity_double[0].header['HISTORY'] = 'r-band stellar mass density map (double resolution)'
mapa_stelar_mass_desity_double.flush()



# PARA EL COLOR (r-i)   (res. doble)************************

counts_double_i=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/i/counts_double_resolution_ri.fits', mode='update')
mu_double_i=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/i/mu_double_resolution_ri.fits', mode='update')
counts_double_r=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/counts_double_resolution_ri.fits', mode='update')


a_lambda=-0.122
b_lambda=1.431
m_sun_i=4.58

log_M_to_L_double=np.arange(len(counts_double_i[0].data[:,0])*1.*len(counts_double_i[0].data[0,:])*1.).reshape(len(counts_double_i[0].data[:,0]),len(counts_double_i[0].data[0,:]))
M_to_L_double=np.arange(len(counts_double_i[0].data[:,0])*1.*len(counts_double_i[0].data[0,:])*1.).reshape(len(counts_double_i[0].data[:,0]),len(counts_double_i[0].data[0,:]))

color_r_i_double=np.arange(len(counts_double_i[0].data[:,0])*1.*len(counts_double_i[0].data[0,:])*1.).reshape(len(counts_double_i[0].data[:,0]),len(counts_double_i[0].data[0,:]))

log_stelar_mass_desity_double=np.arange(len(counts_double_i[0].data[:,0])*1.*len(counts_double_i[0].data[0,:])*1.).reshape(len(counts_double_i[0].data[:,0]),len(counts_double_i[0].data[0,:]))


for x in range(len(counts_double_i[0].data[:,0])):  
  for y in range(len(counts_double_i[0].data[0,:])):
	color_r_i_double[x,y]=-2.5*np.log10(counts_double_r[0].data[x,y]/counts_double_i[0].data[x,y]) # mismo zp para todos los filtros en Stripe82
	log_M_to_L_double[x,y]=(a_lambda+b_lambda*color_r_i_double[x,y])-0.15
	log_stelar_mass_desity_double[x,y]=log_M_to_L_double[x,y]-0.4*(-2.5*np.log10(counts_double_i[0].data[x,y])-m_sun_i)+8.629
		
n=np.arange(len(mu_double_i[0].data[:,0])*1.*len(mu_double_i[0].data[0,:])*1.).reshape(len(mu_double_i[0].data[:,0]),len(mu_double_i[0].data[0,:]))
hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_color_r_i_double.fits')
mapa_color_r_i_double=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_color_r_i_double.fits', mode='update')
mapa_color_r_i_double[0].data[:]=color_r_i_double[:]
mapa_color_r_i_double[0].header=counts_double_i[0].header
mapa_color_r_i_double[0].header['HISTORY'] = '(r-i) color map (double resolution)'
mapa_color_r_i_double.flush()


hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_log_M_to_L_i_double.fits')
mapa_log_M_to_L_double=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_log_M_to_L_i_double.fits', mode='update')
mapa_log_M_to_L_double[0].data[:]=log_M_to_L_double[:]
mapa_log_M_to_L_double[0].header=counts_double_i[0].header
mapa_log_M_to_L_double[0].header['HISTORY'] = 'i-band M/L map (double resolution)'
mapa_log_M_to_L_double.flush()


hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_StellarMassDensity_i_double.fits')
mapa_stelar_mass_desity_double=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_StellarMassDensity_i_double.fits', mode='update')
mapa_stelar_mass_desity_double[0].data[:]=log_stelar_mass_desity_double[:]
mapa_stelar_mass_desity_double[0].header=counts_double_i[0].header
mapa_stelar_mass_desity_double[0].header['HISTORY'] = 'i-band stellar mass density map (double resolution)'
mapa_stelar_mass_desity_double.flush()


##############################################################################################################
##############################################################################################################


##triple RESOLUTION  (mismo proceso pero a resolucion triple)
######################

# PARA EL COLOR (g-r)
counts_triple_r=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/counts_triple_resolution_gr.fits', mode='update')
mu_triple_r=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/mu_triple_resolution_gr.fits', mode='update')
counts_triple_g=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/g/counts_triple_resolution_gr.fits', mode='update')

a_lambda=-0.306
b_lambda=1.097
m_sun_i=4.58

log_M_to_L_triple=np.arange(len(counts_triple_r[0].data[:,0])*1.*len(counts_triple_r[0].data[0,:])*1.).reshape(len(counts_triple_r[0].data[:,0]),len(counts_triple_r[0].data[0,:]))
M_to_L_triple=np.arange(len(counts_triple_r[0].data[:,0])*1.*len(counts_triple_r[0].data[0,:])*1.).reshape(len(counts_triple_r[0].data[:,0]),len(counts_triple_r[0].data[0,:]))

color_g_r_triple=np.arange(len(counts_triple_r[0].data[:,0])*1.*len(counts_triple_r[0].data[0,:])*1.).reshape(len(counts_triple_r[0].data[:,0]),len(counts_triple_r[0].data[0,:]))

log_stelar_mass_desity_triple=np.arange(len(counts_triple_r[0].data[:,0])*1.*len(counts_triple_r[0].data[0,:])*1.).reshape(len(counts_triple_r[0].data[:,0]),len(counts_triple_r[0].data[0,:]))


for x in range(len(counts_triple_r[0].data[:,0])):  
  for y in range(len(counts_triple_r[0].data[0,:])):
	color_g_r_triple[x,y]=-2.5*np.log10(counts_triple_g[0].data[x,y]/counts_triple_r[0].data[x,y]) # mismo zp para todos los filtros en Stripe82
	log_M_to_L_triple[x,y]=(a_lambda+b_lambda*color_g_r_triple[x,y])-0.15
	log_stelar_mass_desity_triple[x,y]=log_M_to_L_triple[x,y]-0.4*(-2.5*np.log10(counts_triple_r[0].data[x,y])-m_sun_r)+8.629
		
n=np.arange(len(mu_triple_r[0].data[:,0])*1.*len(mu_triple_r[0].data[0,:])*1.).reshape(len(mu_triple_r[0].data[:,0]),len(mu_triple_r[0].data[0,:]))
hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_color_g_r_triple.fits')
mapa_color_g_r_triple=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_color_g_r_triple.fits', mode='update')
mapa_color_g_r_triple[0].data[:]=color_g_r_triple[:]
mapa_color_g_r_triple[0].header=counts_triple_r[0].header
mapa_color_g_r_triple[0].header['HISTORY'] = '(g-r) color map (triple resolution)'
mapa_color_g_r_triple.flush()


hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_log_M_to_L_r_triple.fits')
mapa_log_M_to_L_triple=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_log_M_to_L_r_triple.fits', mode='update')
mapa_log_M_to_L_triple[0].data[:]=log_M_to_L_triple[:]
mapa_log_M_to_L_triple[0].header=counts_triple_r[0].header
mapa_log_M_to_L_triple[0].header['HISTORY'] = 'r-band M/L map (triple resolution)'
mapa_log_M_to_L_triple.flush()


hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_StellarMassDensity_r_triple.fits')
mapa_stelar_mass_desity_triple=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_StellarMassDensity_r_triple.fits', mode='update')
mapa_stelar_mass_desity_triple[0].data[:]=log_stelar_mass_desity_triple[:]
mapa_stelar_mass_desity_triple[0].header=counts_triple_r[0].header
mapa_stelar_mass_desity_triple[0].header['HISTORY'] = 'r-band stellar mass density map (triple resolution)'
mapa_stelar_mass_desity_triple.flush()



# PARA EL COLOR (r-i)   (res. triple)************************

counts_triple_i=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/i/counts_triple_resolution_ri.fits', mode='update')
mu_triple_i=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/i/mu_triple_resolution_ri.fits', mode='update')
counts_triple_r=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/r/counts_triple_resolution_ri.fits', mode='update')

a_lambda=-0.122
b_lambda=1.431
m_sun_i=4.58

log_M_to_L_triple=np.arange(len(counts_triple_i[0].data[:,0])*1.*len(counts_triple_i[0].data[0,:])*1.).reshape(len(counts_triple_i[0].data[:,0]),len(counts_triple_i[0].data[0,:]))
M_to_L_triple=np.arange(len(counts_triple_i[0].data[:,0])*1.*len(counts_triple_i[0].data[0,:])*1.).reshape(len(counts_triple_i[0].data[:,0]),len(counts_triple_i[0].data[0,:]))

color_r_i_triple=np.arange(len(counts_triple_i[0].data[:,0])*1.*len(counts_triple_i[0].data[0,:])*1.).reshape(len(counts_triple_i[0].data[:,0]),len(counts_triple_i[0].data[0,:]))

log_stelar_mass_desity_triple=np.arange(len(counts_triple_i[0].data[:,0])*1.*len(counts_triple_i[0].data[0,:])*1.).reshape(len(counts_triple_i[0].data[:,0]),len(counts_triple_i[0].data[0,:]))


for x in range(len(counts_triple_i[0].data[:,0])):  
  for y in range(len(counts_triple_i[0].data[0,:])):
	color_r_i_triple[x,y]=-2.5*np.log10(counts_triple_r[0].data[x,y]/counts_triple_i[0].data[x,y]) # mismo zp para todos los filtros en Stripe82
	log_M_to_L_triple[x,y]=(a_lambda+b_lambda*color_r_i_triple[x,y])-0.15
	log_stelar_mass_desity_triple[x,y]=log_M_to_L_triple[x,y]-0.4*(-2.5*np.log10(counts_triple_i[0].data[x,y])-m_sun_i)+8.629
		
n=np.arange(len(mu_triple_i[0].data[:,0])*1.*len(mu_triple_i[0].data[0,:])*1.).reshape(len(mu_triple_i[0].data[:,0]),len(mu_triple_i[0].data[0,:]))
hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_color_r_i_triple.fits')
mapa_color_r_i_triple=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_color_r_i_triple.fits', mode='update')
mapa_color_r_i_triple[0].data[:]=color_r_i_triple[:]
mapa_color_r_i_triple[0].header=counts_triple_i[0].header
mapa_color_r_i_triple[0].header['HISTORY'] = '(r-i) color map (triple resolution)'
mapa_color_r_i_triple.flush()


hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_log_M_to_L_i_triple.fits')
mapa_log_M_to_L_triple=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_log_M_to_L_i_triple.fits', mode='update')
mapa_log_M_to_L_triple[0].data[:]=log_M_to_L_triple[:]
mapa_log_M_to_L_triple[0].header=counts_triple_i[0].header
mapa_log_M_to_L_triple[0].header['HISTORY'] = 'i-band M/L map (triple resolution)'
mapa_log_M_to_L_triple.flush()


hdu = fits.PrimaryHDU(n)
hdu.writeto('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_StellarMassDensity_i_triple.fits')
mapa_stelar_mass_desity_triple=fits.open('/Users/cris/Tesis/scratch/Stripe82_fitting/UGC01040/mapa_StellarMassDensity_i_triple.fits', mode='update')
mapa_stelar_mass_desity_triple[0].data[:]=log_stelar_mass_desity_triple[:]
mapa_stelar_mass_desity_triple[0].header=counts_triple_i[0].header
mapa_stelar_mass_desity_triple[0].header['HISTORY'] = 'i-band stellar mass density map (triple resolution)'
mapa_stelar_mass_desity_triple.flush()


##############################################################################################################
##############################################################################################################

"""
# 2.  total area of the galaxy -> Total Mass

Dpc=DMpc*1e6
semi_mayor_pc=math.tan(((semi_mayor*scale)/(60.*60.))*(math.pi/180.))*Dpc
semi_menor_pc=math.tan(((semi_menor*scale)/(60.*60.))*(math.pi/180.))*Dpc

#galaxy_sup=math.pi*semi_mayor_pc*semi_menor_pc   # tengo que ver como hacer esto...solo tengo que coger los pixeles pertenecientes a la galaxia, quiza segun algun criterio en la magnitud.

#total_mass_Msun=np.nansum(10**log_stelar_mass_desity)*galaxy_sup


# 3.  comparison with dynamical mass

#dyn_mass=2.326e5*(np.log10(log_d25)/2.)*v_rot**2   ***logd25 debo cambiarlo a las uds que son (me lo dan en 0.1arcsec o algo asi...)
"""

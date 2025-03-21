import numpy as np
import os 
from astropy.io import fits
import matplotlib.pyplot as plt
plt.ion()
import math
from astropy.stats.funcs import biweight_location
from numpy import mean, sqrt, square, arange
from matplotlib import rc, rcParams
from matplotlib import gridspec



modSSP=np.loadtxt('/Users/cris/Tesis/Mod_EMILES_SSP/mag_tmpRHvWmg_OldPop.mag')
NUV_old_mag = modSSP[:,2]
r_old_mag = modSSP[:,3]
spitz36m_old_mag = modSSP[:,4]


percent_old =1.
percent_young = 1. - percent_old

#age = 0.05 #Gyr
#age = 0.1 #Gyr
age = 0.2 #Gyr
#age = 0.3 #Gyr

#age = 0.4 #Gyr
#age = 0.5 #Gyr
#age = 0.6 #Gyr


# Paso a Jy toda la tabla de poblacion vieja y multiplico por 0.995:

NUV_old_Jy = 10.**((-2./5.)*(NUV_old_mag - 8.9)) * percent_old
r_old_Jy = 10.**((-2./5.)*(r_old_mag - 8.9)) * percent_old
spitz36m_old_Jy = 10.**((-2./5.)*(spitz36m_old_mag - 8.9)) * percent_old



# Poblacion Joven: metalicidad solar, edad 0.1 Gyr -- paso a Jy y multiplico por 0.005
if age == 0.05:
	NUV_01Gyr_Jy = 10.**((-2./5.)*(3.24469 - 8.9)) * percent_young
	r_01Gyr_Jy = 10.**((-2./5.)*(2.36608 - 8.9)) * percent_young
	spitz36m_01Gyr_Jy = 10.**((-2./5.)*(-99. - 8.9)) * percent_young

if age == 0.1:
	NUV_01Gyr_Jy = 10.**((-2./5.)*(4.14753 - 8.9)) * percent_young
	r_01Gyr_Jy = 10.**((-2./5.)*(2.89035 - 8.9)) * percent_young
	spitz36m_01Gyr_Jy = 10.**((-2./5.)*(3.49993 - 8.9)) * percent_young

if age == 0.2:
	NUV_01Gyr_Jy = 10.**((-2./5.)*(5.08906 - 8.9)) * percent_young
	r_01Gyr_Jy = 10.**((-2./5.)*(3.35353 - 8.9)) * percent_young
	spitz36m_01Gyr_Jy = 10.**((-2./5.)*(4.12326 - 8.9)) * percent_young

if age == 0.3:
	NUV_01Gyr_Jy = 10.**((-2./5.)*(5.81193 - 8.9)) * percent_young
	r_01Gyr_Jy = 10.**((-2./5.)*(3.62002 - 8.9)) * percent_young
	spitz36m_01Gyr_Jy = 10.**((-2./5.)*(4.34592 - 8.9)) * percent_young

if age == 0.4:
	NUV_01Gyr_Jy = 10.**((-2./5.)*(6.25683 - 8.9)) * percent_young
	r_01Gyr_Jy = 10.**((-2./5.)*(3.83112 - 8.9)) * percent_young
	spitz36m_01Gyr_Jy = 10.**((-2./5.)*(4.47385 - 8.9)) * percent_young

if age == 0.5:
	NUV_01Gyr_Jy = 10.**((-2./5.)*(6.81075 - 8.9)) * percent_young
	r_01Gyr_Jy = 10.**((-2./5.)*(3.99947 - 8.9)) * percent_young
	spitz36m_01Gyr_Jy = 10.**((-2./5.)*(4.58480 - 8.9)) * percent_young

if age == 0.6:
	NUV_01Gyr_Jy = 10.**((-2./5.)*(7.30911 - 8.9)) * percent_young
	r_01Gyr_Jy = 10.**((-2./5.)*(4.13603 - 8.9)) * percent_young
	spitz36m_01Gyr_Jy = 10.**((-2./5.)*(4.66622 - 8.9)) * percent_young


# Sumo ambas contribuciones
NUV_Jy = NUV_old_Jy + NUV_01Gyr_Jy
r_Jy = r_old_Jy + r_01Gyr_Jy 
spitz36m_Jy = spitz36m_old_Jy + spitz36m_01Gyr_Jy

# Paso a AB mag

NUV_AB = (-5./2.)*np.log10(NUV_Jy/3631.)
r_AB = (-5./2.)*np.log10(r_Jy/3631.)
spitz36m_AB = (-5./2.)*np.log10(spitz36m_Jy/3631.)


# Guardo los datos
file_mags=np.loadtxt('/Users/cris/Tesis/Mod_EMILES_SSP/mag_tmpRHvWmg.mag')
file_mags[:,2] = NUV_AB
file_mags[:,3] = r_AB
file_mags[:,4] = spitz36m_AB

np.savetxt('/Users/cris/Tesis/Mod_EMILES_SSP/mag_tmpRHvWmg.mag', file_mags)






import numpy as np
import os 
from astropy.table import Table

import matplotlib.pyplot as plt
plt.ion()
plt.close('all')

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14

#plt.rc('text', usetex=True)

# =============================================================================
# Parametros iniciales del grupo/imagen/IGL/aperturas
# =============================================================================
groupID = '400138'
members = ['IGL']

gal_lowrlim = 10 
gal_uprlim = 60 # pix; Rmax up to where obtain the colors of the galaxies
IGL_lowrlim = 90 # pix
IGL_uprlim = 210 # pix

# =============================================================================
# Paths
# =============================================================================
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'



# =============================================================================
# Creamos directorios para guardar los distintos productos del programa
# =============================================================================

try:
    # Create target Directory
    os.mkdir(results_path+'IGL/StellPop')
except:
    pass

# =============================================================================
# Fichero de parametros
# =============================================================================
# Copiamos el fichero de parametros para asegurarnos de que es el mismo

cmd = 'cp '+results_path+'params.py .'

os.system(cmd)
import params

# Group params
bcgID = params.bcgID
z = params.z


# Instrument params
zp = params.zp
pixscale = params.pixelscale

# Extinction and K-corrs
k_corrs = params.k_corrs
extinctions = params.extinctions
dimming = 10*np.log10(1+z) 

   
os.system('rm params.py')


# =============================================================================
#Files with the color profiles
# =============================================================================

colors = ['g-i', 'g-r','r-i']

files = {'g-i':[results_path+'IGL/color/profiles/IGL_prof_g-i.fits', results_path+'IGL/color/profiles/group_prof_g-i.fits'],
         'g-r':[results_path+'IGL/color/profiles/IGL_prof_g-r.fits', results_path+'IGL/color/profiles/group_prof_g-r.fits'], 
         'r-i':[results_path+'IGL/color/profiles/IGL_prof_r-i.fits', results_path+'IGL/color/profiles/group_prof_r-i.fits']}

   

# =============================================================================
# Load the Photometric Predictions based on E-MILES SEDs
# =============================================================================
 
# MODELO DE POBLACIONES ESTELARES VAZDEKIS
# The appropriate references for these predictions are Vazdekis et al. (2012), Ricciardelli et al. (2012) and Vazdekis et al. (2016). 
# http://research.iac.es/proyecto/miles/pages/photometric-predictions-based-on-e-miles-seds.php

modVaz=np.loadtxt('/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/data/EMILES-SEDs/sdss_ch_iPp0.00.MAG')


metallicity = modVaz[:,1]
age = modVaz[:,2]
g_SDSS = modVaz[:,5]
r_SDSS = modVaz[:,6]
i_SDSS = modVaz[:,7]
ML_g_SDSS = modVaz[:,10]
ML_r_SDSS = modVaz[:,11]
ML_i_SDSS = modVaz[:,12]
 
metall_23 = modVaz[0:50,1]
metall_17 = modVaz[50:100,1]
metall_13 = modVaz[100:150,1]
metall_07 = modVaz[150:200,1]
metall_04 = modVaz[200:250,1]
metall_00 = modVaz[250:300,1]
metall_02 = modVaz[300:350,1]

age_23 = modVaz[0:50,2]
age_17 = modVaz[50:100,2]
age_13 = modVaz[100:150,2]
age_07 = modVaz[150:200,2]
age_04 = modVaz[200:250,2]
age_00 = modVaz[250:300,2]
age_02 = modVaz[300:350,2]

g_r_SDSS = g_SDSS - r_SDSS
g_i_SDSS = g_SDSS - i_SDSS
r_i_SDSS = r_SDSS - i_SDSS

g_r_SDSS_23 = g_r_SDSS[0:50]
g_r_SDSS_17 = g_r_SDSS[50:100]
g_r_SDSS_13 = g_r_SDSS[100:150]
g_r_SDSS_07 = g_r_SDSS[150:200]
g_r_SDSS_04 = g_r_SDSS[200:250]
g_r_SDSS_00 = g_r_SDSS[250:300]
g_r_SDSS_02 = g_r_SDSS[300:350]


g_i_SDSS_23 = g_i_SDSS[0:50]
g_i_SDSS_17 = g_i_SDSS[50:100]
g_i_SDSS_13 = g_i_SDSS[100:150]
g_i_SDSS_07 = g_i_SDSS[150:200]
g_i_SDSS_04 = g_i_SDSS[200:250]
g_i_SDSS_00 = g_i_SDSS[250:300]
g_i_SDSS_02 = g_i_SDSS[300:350]

r_i_SDSS_23 = r_i_SDSS[0:50]
r_i_SDSS_17 = r_i_SDSS[50:100]
r_i_SDSS_13 = r_i_SDSS[100:150]
r_i_SDSS_07 = r_i_SDSS[150:200]
r_i_SDSS_04 = r_i_SDSS[200:250]
r_i_SDSS_00 = r_i_SDSS[250:300]
r_i_SDSS_02 = r_i_SDSS[300:350]



# =============================================================================
# COLOR - COLOR PLOTS g-r vs r-i   +  Stellar Populations VAZDK
# =============================================================================


# Leemos los perfiles que necesitamos: 
# g-r
g_r_IGL = files['g-r'][0]

r = Table.read(g_r_IGL)['r'] # La posicion radial es la misma para todos los perfiles
color_g_r = Table.read(g_r_IGL)['color']
err_g_r = Table.read(g_r_IGL)['err']


g_r_IGLGal = files['g-r'][1]

r_gal = Table.read(g_r_IGLGal)['r'] # La posicion radial es la misma para todos los perfiles
color_g_r_gal = Table.read(g_r_IGLGal)['color']
err_g_r_gal = Table.read(g_r_IGLGal)['err']


# r-i
r_i_IGL = files['r-i'][0]

r = Table.read(r_i_IGL)['r'] # La posicion radial es la misma para todos los perfiles
color_r_i = Table.read(r_i_IGL)['color']
err_r_i = Table.read(r_i_IGL)['err']

r_i_IGLGal = files['r-i'][1]

r_gal = Table.read(r_i_IGLGal)['r'] # La posicion radial es la misma para todos los perfiles
color_r_i_gal = Table.read(r_i_IGLGal)['color']
err_r_i_gal = Table.read(r_i_IGLGal)['err']


# *****************************************************************************

# Plot
plt.figure(3,figsize=(11,7.5))
plt.grid(3)

cmap = plt.get_cmap('coolwarm')
colors = cmap(np.linspace(0, 0.9, 7))

# Ages points
ages_in = [24, 32, 38, 44] # 0, 4, 12, 18, 24, 32, 38, 44
ages_label = ['1 Gyr', '2.5 Gyr', '5 Gyr', '10 Gyr'] # '0.06 Gyr', '0.1 Gyr','0.25 Gyr', '0.5 Gyr', '1 Gyr', '2.5 Gyr', '5 Gyr', '10 Gyr'
ages_ms = np.linspace(4, 25, len(ages_in))


# Metallicity tracks
#plt.plot(g_r_SDSS_23[ages_in[0]:ages_in[-1]+1],r_i_SDSS_23[ages_in[0]:ages_in[-1]+1], lw=3, color=colors[0],  label='-2.3 [Fe/H]', zorder=1)
#plt.plot(g_r_SDSS_17[ages_in[0]:ages_in[-1]+1],r_i_SDSS_17[ages_in[0]:ages_in[-1]+1], lw=3, color=colors[1], label='-1.7', zorder=1)
plt.plot(g_r_SDSS_13[ages_in[0]:ages_in[-1]+1],r_i_SDSS_13[ages_in[0]:ages_in[-1]+1], lw=3, color=colors[2], label='-1.3 [Fe/H]', zorder=1)
plt.plot(g_r_SDSS_07[ages_in[0]:ages_in[-1]+1],r_i_SDSS_07[ages_in[0]:ages_in[-1]+1], lw=3, color=colors[3], label='-0.7', zorder=1)
plt.plot(g_r_SDSS_04[ages_in[0]:ages_in[-1]+1],r_i_SDSS_04[ages_in[0]:ages_in[-1]+1], lw=3, color=colors[4], label='-0.4', zorder=1)
plt.plot(g_r_SDSS_00[ages_in[0]:ages_in[-1]+1],r_i_SDSS_00[ages_in[0]:ages_in[-1]+1], lw=3, color=colors[5], label='0.0', zorder=1)
plt.plot(g_r_SDSS_02[ages_in[0]:ages_in[-1]+1],r_i_SDSS_02[ages_in[0]:ages_in[-1]+1], lw=3, color=colors[6], label='0.2', zorder=1)


for ind, ind_table in enumerate(ages_in):
    
    # Puntos para las edades
    #plt.plot(g_r_SDSS_23[ind_table],r_i_SDSS_23[ind_table],'o', markersize=ages_ms[ind], color=colors[0], alpha=0.6, zorder=2,  label=ages_label[ind])
    #plt.plot(g_r_SDSS_17[ind_table],r_i_SDSS_17[ind_table],'o', markersize=ages_ms[ind], color=colors[1], alpha=0.6, zorder=2)
    plt.plot(g_r_SDSS_13[ind_table],r_i_SDSS_13[ind_table],'o', markersize=ages_ms[ind], color=colors[2], alpha=0.6, zorder=2, label=ages_label[ind])
    plt.plot(g_r_SDSS_07[ind_table],r_i_SDSS_07[ind_table],'o', markersize=ages_ms[ind], color=colors[3], alpha=0.6, zorder=2)
    plt.plot(g_r_SDSS_04[ind_table],r_i_SDSS_04[ind_table],'o', markersize=ages_ms[ind], color=colors[4], alpha=0.6, zorder=2)
    plt.plot(g_r_SDSS_00[ind_table],r_i_SDSS_00[ind_table],'o', markersize=ages_ms[ind], color=colors[5], alpha=0.6, zorder=2)
    plt.plot(g_r_SDSS_02[ind_table],r_i_SDSS_02[ind_table],'o', markersize=ages_ms[ind], color=colors[6], alpha=0.6, zorder=2)



'''
# My data colors
last_point = -4
data_ms = np.linspace(3, 30, len(color_g_r[:last_point]))

# points with increasing size
for ii in np.arange(7, len(color_g_r[:last_point])):
    plt.errorbar(color_g_r[ii], color_r_i[ii], xerr=err_g_r[ii], yerr=err_r_i[ii], marker='*',linestyle='None', markersize=data_ms[ii], color='k', alpha = 0.9, zorder=3)
    #plt.errorbar(color_g_r_gal[ii], color_r_i_gal[ii], xerr=err_g_r_gal[ii], yerr=err_r_i_gal[ii], marker='^',linestyle='None', markersize=data_ms[ii], color='y', alpha = 0.9, zorder=3)

'''
# =============================================================================
# Plot mean values for the IGL AND galaxies from my data
# =============================================================================


IGL_region = (r>IGL_lowrlim) & (r<IGL_uprlim)
gal_region = (r>gal_lowrlim) & (r<gal_uprlim)

mean_IGLcolor_g_r = np.nanmean(color_g_r[IGL_region])
mean_IGLcolor_r_i= np.nanmean(color_r_i[IGL_region])

mean_GALcolor_g_r = np.nanmean(color_g_r_gal[gal_region])
mean_GALcolor_r_i= np.nanmean(color_r_i_gal[gal_region])


# THE ERROR OF THE MEAN IS THE STD OF THE ERRORS DISTRIBUTION 
mean_IGLerr_g_r = np.nanstd(err_g_r[IGL_region])
mean_IGLerr_r_i = np.nanstd(err_r_i[IGL_region])

mean_GALerr_g_r = np.nanstd(err_g_r_gal[gal_region])
mean_GALerr_r_i = np.nanstd(err_r_i_gal[gal_region])

# Plot
plt.errorbar(mean_IGLcolor_g_r, mean_IGLcolor_r_i, xerr=mean_IGLerr_g_r, yerr=mean_IGLerr_r_i, marker='X', markersize=10, linestyle='None', color='C8', zorder=3, alpha = 0.9, label='IGL')
plt.errorbar(mean_GALcolor_g_r, mean_GALcolor_r_i, xerr=mean_GALerr_g_r, yerr=mean_GALerr_r_i, marker='X',linestyle='None', markersize=10, color='C7', zorder=3, alpha = 0.9, label='Group')





# =============================================================================
# Plot mean values for the galaxies from CFHT data -- SANITY CHECK!
# =============================================================================
CFHT_g_mag = np.array([18.77734, 19.51939, 18.861282])
CFHT_r_mag = np.array([17.64407, 18.387314, 17.73466])
CFHT_i_mag = np.array([17.094957, 17.837996, 17.185776])

CFHT_g_magerr = np.mean([6.3821004E-4, 8.563916E-4, 7.7374425E-4])
CFHT_r_magerr = np.mean([4.0916473E-4, 5.750851E-4, 4.8224913E-4])
CFHT_i_magerr = np.mean([3.554823E-4, 5.0398847E-4, 4.4499396E-4])

#CFHT_g_ext = 0.0894389
#CFHT_r_ext = 0.0648686
#CFHT_i_ext = 0.0491879


CFHT_g = CFHT_g_mag - k_corrs['g']
CFHT_r = CFHT_r_mag - k_corrs['r']
CFHT_i = CFHT_i_mag - k_corrs['i']

CFHT_g_r = CFHT_g-CFHT_r
CFHT_r_i = CFHT_r-CFHT_i

CFHT_g_r_err = np.sqrt(np.square(CFHT_g_magerr) + np.square(CFHT_r_magerr))
CFHT_r_i_err = np.sqrt(np.square(CFHT_r_magerr) + np.square(CFHT_i_magerr))


#plt.errorbar(np.mean(CFHT_g_r), np.mean(CFHT_r_i), xerr=CFHT_g_r_err, yerr=CFHT_r_i_err, marker='D', markersize=3, linestyle='None', color='darkgreen', zorder=3, alpha = 0.9, label='Group CFHT')

SDSS_g_mag = np.array([19.986774, 20.89709, 20.33445])
SDSS_r_mag = np.array([18.895025, 19.734497, 19.162361])
SDSS_i_mag = np.array([18.3107, 19.111395, 18.590382])

SDSS_g_magerr = np.std([0.032168753, 0.08190984, 0.0376867])
SDSS_r_magerr = np.std([0.028846331, 0.07544215, 0.031611435])
SDSS_i_magerr = np.std([0.028017152, 0.07069266, 0.031041129])

#SDSS_g_ext = 0.09149324
#SDSS_r_ext = 0.06635853
#SDSS_i_ext = 0.050317664

SDSS_g = SDSS_g_mag - k_corrs['g']
SDSS_r = SDSS_r_mag - k_corrs['r']
SDSS_i = SDSS_i_mag - k_corrs['i']

SDSS_g_r = SDSS_g-SDSS_r
SDSS_r_i = SDSS_r-SDSS_i

SDSS_g_r_err = np.sqrt(np.square(SDSS_g_magerr) + np.square(SDSS_r_magerr))
SDSS_r_i_err = np.sqrt(np.square(SDSS_r_magerr) + np.square(SDSS_i_magerr))



#plt.errorbar(np.mean(SDSS_g_r), np.mean(SDSS_r_i), xerr=SDSS_g_r_err, yerr=SDSS_r_i_err, marker='s', markersize=8, linestyle='None', color='limegreen', zorder=3, alpha = 0.9, label='Group SDSS')




plt.ylabel("$r-i$")
plt.xlabel("$g-r$")
#plt.xlim(-0.1,0.8)
#plt.ylim(-0.1,0.4)

plt.legend(loc=2, fontsize=12, numpoints=1, ncol=3)	

plt.savefig(results_path+'IGL/StellPop/Color_Color.pdf')























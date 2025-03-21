import numpy as np
import os 
from astropy.io import fits
import matplotlib.pyplot as plt
plt.ion()
from numpy import mean, sqrt, square, arange
from matplotlib import rc, rcParams
plt.close('all')

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14

plt.rc('text', usetex=True)



g_r_Thin_data = np.loadtxt('/Users/cris/Tesis/scratch/Stripe82_fitting/'+galaxy+'/Color_DataProf_ThinDisk_g-r.txt')		
g_i_Thin_data = np.loadtxt('/Users/cris/Tesis/scratch/Stripe82_fitting/'+galaxy+'/Color_DataProf_ThinDisk_g-i.txt')		
r_i_Thin_data = np.loadtxt('/Users/cris/Tesis/scratch/Stripe82_fitting/'+galaxy+'/Color_DataProf_ThinDisk_r-i.txt')		


# MODELO DE POBLACIONES ESTELARES VAZDEKIS
# The appropriate references for these predictions are Vazdekis et al. (2012), Ricciardelli et al. (2012) and Vazdekis et al. (2016). 
# http://research.iac.es/proyecto/miles/pages/photometric-predictions-based-on-e-miles-seds.php

modVaz=np.loadtxt('/Users/cris/Tesis/Mod_EMILES_SSP/sdss_ch_iPp0.00.MAG')
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





# ***************************************************************
# EMPIEZO CON LOS PLOTS PARA CADA GALAXIA
# ***************************************************************



# ***************************************************************
# COLOR - COLOR PLOTS g-r vs r-i   + POB ESTELARES VAZDK
# ***************************************************************


size_td = [1.+n**3.5 for n in  list(range(len(g_r_Thin_deconv[0])))]
size_Td = [1.+n**3.5 for n in  list(range(len(g_r_Thick_deconv[0])))]



#*****************************************************************
# THIN DISC
#*****************************************************************

plt.figure(3,figsize=(11,7.5))
plt.grid(3)

cmap = plt.get_cmap('coolwarm')
colors = cmap(np.linspace(0, 1.0, 7))


plt.plot(g_r_SDSS_23,r_i_SDSS_23, 'o-', lw=3, color=colors[0], mec=colors[0], label='-2.3 [Fe/H]', zorder=1)
plt.plot(g_r_SDSS_17,r_i_SDSS_17, 'o-', lw=3, color=colors[1], mec=colors[1], label='		 -1.7', zorder=1)
plt.plot(g_r_SDSS_13,r_i_SDSS_13, 'o-', lw=3, color=colors[2], mec=colors[2], label='       -1.3', zorder=1)
plt.plot(g_r_SDSS_07,r_i_SDSS_07 , 'o-', lw=3, color=colors[3], mec=colors[3], label='       -0.7', zorder=1)
plt.plot(g_r_SDSS_04,r_i_SDSS_04 , 'o-', lw=3, color=colors[4], mec=colors[4], label='       -0.4', zorder=1)
plt.plot(g_r_SDSS_00,r_i_SDSS_00 , 'o-', lw=3, color=colors[5], mec=colors[5], label='        0.0', zorder=1)
plt.plot(g_r_SDSS_02,r_i_SDSS_02 , 'o-', lw=3, color=colors[6], mec=colors[6], label='        0.2', zorder=1)

plt.scatter(0, 0, color='w', alpha=0.1, label=' ')


# Puntos para las edades
plt.plot(g_r_SDSS_23[18],r_i_SDSS_23[18],'s', markersize=15, color='m', mec='m', zorder=2,  label='0.5 Gyr')
plt.plot(g_r_SDSS_17[18],r_i_SDSS_17[18],'s', markersize=15, color='m', mec='m', zorder=2)
plt.plot(g_r_SDSS_13[18],r_i_SDSS_13[18],'s', markersize=15, color='m', mec='m', zorder=2)
plt.plot(g_r_SDSS_07[18],r_i_SDSS_07[18],'s', markersize=15, color='m', mec='m', zorder=2)
plt.plot(g_r_SDSS_04[18],r_i_SDSS_04[18],'s', markersize=15, color='m', mec='m', zorder=2)
plt.plot(g_r_SDSS_00[18],r_i_SDSS_00[18],'s', markersize=15, color='m', mec='m', zorder=2)
plt.plot(g_r_SDSS_02[18],r_i_SDSS_02[18],'s', markersize=15, color='m', mec='m', zorder=2)

plt.plot(g_r_SDSS_23[24],r_i_SDSS_23[24],'*', markersize=20, color='gold', mec='gold', zorder=2, label='1 Gyr')
plt.plot(g_r_SDSS_17[24],r_i_SDSS_17[24],'*', markersize=20, color='gold', mec='gold', zorder=2)
plt.plot(g_r_SDSS_13[24],r_i_SDSS_13[24],'*', markersize=20, color='gold', mec='gold', zorder=2)
plt.plot(g_r_SDSS_07[24],r_i_SDSS_07[24],'*', markersize=20, color='gold', mec='gold', zorder=2)
plt.plot(g_r_SDSS_04[24],r_i_SDSS_04[24],'*', markersize=20, color='gold', mec='gold', zorder=2)
plt.plot(g_r_SDSS_00[24],r_i_SDSS_00[24],'*', markersize=20, color='gold', mec='gold', zorder=2)
plt.plot(g_r_SDSS_02[24],r_i_SDSS_02[24],'*', markersize=20, color='gold', mec='gold', zorder=2)

plt.plot(g_r_SDSS_23[32],r_i_SDSS_23[32],'^', markersize=15, color='b', mec='b', zorder=2, label='2.5 Gyr')
plt.plot(g_r_SDSS_17[32],r_i_SDSS_17[32],'^', markersize=15, color='b', mec='b', zorder=2)
plt.plot(g_r_SDSS_13[32],r_i_SDSS_13[32],'^', markersize=15, color='b', mec='b', zorder=2)
plt.plot(g_r_SDSS_07[32],r_i_SDSS_07[32],'^', markersize=15, color='b', mec='b', zorder=2)
plt.plot(g_r_SDSS_04[32],r_i_SDSS_04[32],'^', markersize=15, color='b', mec='b', zorder=2)
plt.plot(g_r_SDSS_00[32],r_i_SDSS_00[32],'^', markersize=15, color='b', mec='b', zorder=2)
plt.plot(g_r_SDSS_02[32],r_i_SDSS_02[32],'^', markersize=15, color='b', mec='b', zorder=2)

plt.plot(g_r_SDSS_23[38],r_i_SDSS_23[38],'h', markersize=15, color='mediumseagreen', mec='mediumseagreen', zorder=2, label='5 Gyr')
plt.plot(g_r_SDSS_17[38],r_i_SDSS_17[38],'h', markersize=15, color='mediumseagreen', mec='mediumseagreen', zorder=2)
plt.plot(g_r_SDSS_13[38],r_i_SDSS_13[38],'h', markersize=15, color='mediumseagreen', mec='mediumseagreen', zorder=2)
plt.plot(g_r_SDSS_07[38],r_i_SDSS_07[38],'h', markersize=15, color='mediumseagreen', mec='mediumseagreen', zorder=2)
plt.plot(g_r_SDSS_04[38],r_i_SDSS_04[38],'h', markersize=15, color='mediumseagreen', mec='mediumseagreen', zorder=2)
plt.plot(g_r_SDSS_00[38],r_i_SDSS_00[38],'h', markersize=15, color='mediumseagreen', mec='mediumseagreen', zorder=2)
plt.plot(g_r_SDSS_02[38],r_i_SDSS_02[38],'h', markersize=15, color='mediumseagreen', mec='mediumseagreen', zorder=2)

plt.plot(g_r_SDSS_23[44],r_i_SDSS_23[44],'v', markersize=15, color='darkorange', mec='darkorange', zorder=2, label='10 Gyr')
plt.plot(g_r_SDSS_17[44],r_i_SDSS_17[44],'v', markersize=15, color='darkorange', mec='darkorange', zorder=2)
plt.plot(g_r_SDSS_13[44],r_i_SDSS_13[44],'v', markersize=15, color='darkorange', mec='darkorange', zorder=2)
plt.plot(g_r_SDSS_07[44],r_i_SDSS_07[44],'v', markersize=15, color='darkorange', mec='darkorange', zorder=2)
plt.plot(g_r_SDSS_04[44],r_i_SDSS_04[44],'v', markersize=15, color='darkorange', mec='darkorange', zorder=2)
plt.plot(g_r_SDSS_00[44],r_i_SDSS_00[44],'v', markersize=15, color='darkorange', mec='darkorange', zorder=2)
plt.plot(g_r_SDSS_02[44],r_i_SDSS_02[44],'v', markersize=15, color='darkorange', mec='darkorange', zorder=2)


plt.errorbar(g_r_Thin_deconv[0],r_i_Thin_deconv[0],xerr=[g_r_Thin_deconv[1],g_r_Thin_deconv[2]], yerr=[r_i_Thin_deconv[1],r_i_Thin_deconv[2]], marker='.',linestyle='None', color='k', zorder=3)
plt.scatter(g_r_Thin_deconv[0],r_i_Thin_deconv[0], s=size_td, edgecolor='None', color='k', zorder=3)
plt.ylabel("$r-i$")
plt.title('%s Thin disc [decon. model]' %(galaxy))
plt.xlabel("$g-r$")
plt.xlim(0.2,0.8)
plt.ylim(0.05,0.35)

plt.legend(loc=4, fontsize=16, numpoints=1, ncol=2)	






plt.savefig('/Users/cris/Tesis/scratch/Stripe82_fitting/'+galaxy+'/Color_Color_Deconv_ThinDisk.png')




'''
# Puntos para las edades
plt.plot(g_r_SDSS_00[38],r_i_SDSS_00[38],'s', markersize=15, color='mediumseagreen',mec='mediumseagreen', label='5 Gyr')
plt.plot(g_r_SDSS_04[38],r_i_SDSS_04[38],'s', markersize=15, color='mediumseagreen', mec='mediumseagreen')
plt.plot(g_r_SDSS_07[38],r_i_SDSS_07[38],'s', markersize=15, color='mediumseagreen', mec='mediumseagreen')

plt.plot(g_r_SDSS_00[41],r_i_SDSS_00[41],'*', markersize=20, color='m', mec='m', label='7 Gyr')
plt.plot(g_r_SDSS_04[41],r_i_SDSS_04[41],'*', markersize=20, color='m', mec='m')
plt.plot(g_r_SDSS_07[41],r_i_SDSS_07[41],'*', markersize=20, color='m', mec='m')

plt.plot(g_r_SDSS_00[44],r_i_SDSS_00[44],'^', markersize=15, color='b', mec='b', label='10 Gyr')
plt.plot(g_r_SDSS_04[44],r_i_SDSS_04[44],'^', markersize=15, color='b', mec='b')
plt.plot(g_r_SDSS_07[44],r_i_SDSS_07[44],'^', markersize=15, color='b', mec='b')


plt.text(g_r_SDSS_00[38],r_i_SDSS_00[38]-0.01,'5 Gyr',color='mediumseagreen', fontweight='heavy')
plt.text(g_r_SDSS_04[38],r_i_SDSS_04[38]+0.01,'5 Gyr',color='mediumseagreen', weight='heavy')
plt.text(g_r_SDSS_07[38],r_i_SDSS_07[38]-0.01,'5 Gyr', color='mediumseagreen', weight='heavy')

plt.text(g_r_SDSS_00[41],r_i_SDSS_00[41]-0.01,'7 Gyr', color='m', weight='heavy')
plt.text(g_r_SDSS_04[41],r_i_SDSS_04[41]+0.01,'7 Gyr', color='m', weight='heavy')
plt.text(g_r_SDSS_07[41],r_i_SDSS_07[41]-0.01,'7 Gyr', color='m', weight='heavy')

plt.text(g_r_SDSS_00[44]-0.015,r_i_SDSS_00[44]-0.01,'10 Gyr', color='b', weight='heavy')
plt.text(g_r_SDSS_04[44],r_i_SDSS_04[44]+0.01,'10 Gyr', color='b', weight='heavy')
plt.text(g_r_SDSS_07[44],r_i_SDSS_07[44]-0.008,'10 Gyr', color='b', weight='heavy')
'''


#plt.xlim(0.5,0.82)
#plt.ylim(0.2,0.42)





'''
# ***************************************************************
# COLOR - COLOR PLOTS g-i vs g-r
# ***************************************************************

size_td = [1.+n**3.5 for n in  list(reversed(range(len(g_r_Thin_deconv[0]))))]
size_Td = [1.+n**3.5 for n in  list(reversed(range(len(g_r_Thick_deconv[0]))))]

plt.figure(1,figsize=(11,7.5))
plt.grid(1)
plt.errorbar(g_i_Thin_deconv[0],g_r_Thin_deconv[0],xerr=[g_i_Thin_deconv[1],g_i_Thin_deconv[2]], yerr=[g_r_Thin_deconv[1],g_r_Thin_deconv[2]], marker='.',linestyle='None')
plt.scatter(g_i_Thin_deconv[0],g_r_Thin_deconv[0], s=size_td, edgecolor='None')
plt.ylabel("$g-r$")
plt.title('%s Thin disc colour-colour [decon. model]' %(galaxy))
plt.xlabel("$g-i$")
#plt.xlim(0.5,1.5)
#plt.ylim(0.5,1.5)
#plt.savefig('/Users/cris/Tesis/scratch/Stripe82_fitting/'+galaxy+'/Color_Color_DeconvProf_ThinDisk.eps')

plt.figure(2,figsize=(11,7.5))
plt.grid(2)
plt.errorbar(g_i_Thick_deconv[0],g_r_Thick_deconv[0],xerr=[g_i_Thick_deconv[1],g_i_Thick_deconv[2]], yerr=[g_r_Thick_deconv[1],g_r_Thick_deconv[2]], marker='.',linestyle='None')
plt.scatter(g_i_Thick_deconv[0],g_r_Thick_deconv[0], s=size_Td, edgecolor='None')
plt.ylabel("$g-r$")
plt.title('%s Thick disc colour-colour [decon. model]' %(galaxy))
plt.xlabel("$g-i$")
#plt.xlim(0.5,1.5)
#plt.ylim(0.5,1.5)
#plt.savefig('/Users/cris/Tesis/scratch/Stripe82_fitting/'+galaxy+'/Color_Color_DeconvProf_ThickDisk.eps')
'''















'''
# ***************************************************************
# COLOR - COLOR "FAKE" PLOTS 
# ***************************************************************

dataThin1 = np.genfromtxt('/Users/cris/Tesis/scratch/Stripe82_fitting/'+galaxy+'/g/Color_DataProf_thinkDisk_g.txt')
x_data = dataThin1[0,:]


plt.figure(5)
plt.grid(5)
plt.plot(x_data,g_r_Thin_deconv[0]-g_r_Thick_deconv[0], marker='.',linestyle='None')
plt.ylabel("g-r")
plt.title('%s Thin-thick disk color g-r [decon. model]' %(galaxy))
plt.xlabel("R [arcsec]")

plt.figure(6)
plt.grid(6)
plt.plot(x_data,g_i_Thin_deconv[0]-g_i_Thick_deconv[0], marker='.',linestyle='None')
plt.ylabel("g-i")
plt.title('%s Thin-thick disk color g-i[decon. model]' %(galaxy))
plt.xlabel("R [arcsec]")

'''
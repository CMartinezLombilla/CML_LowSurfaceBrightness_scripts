import numpy as np
import numpy.ma as ma
import os
import matplotlib.pyplot as plt
plt.close('all')
plt.ioff()

plt.rcParams["font.family"] = "serif"
plt.rc('text', usetex=False)

from astropy.io import fits



# =============================================================================
# importamos CML
# =============================================================================

import sys
sys.path.append('/Users/felipe/Dropbox/CML/') # path del CML Feli
sys.path.append('/Users/z3530031/Dropbox/CML/') # path del CML Cris

sys.path.append('/Users/felipe/Dropbox/CML/WHT/') # path del fichero de parametros Feli
sys.path.append('/Users/z3530031/Dropbox/CML/WHT/') # path del fichero de parametros Cris

from CML import EllipticalProfile
from CML import LogNorm
from CML import Pix2Kpc
from CML import Kpc2Pix

# Division entera en "base"
def baseround(x, base=10):
    return base * round(x/base)

#==============================================================================
# Parametros iniciales
#==============================================================================
GroupID = '400138'
member = 'BCG-sky'

# Fijar el centro de las elipses
fix_center = True


# Banda a utilizar en el analisis principal
band = 'i'   

# =============================================================================
# Paths
# =============================================================================
masks_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+GroupID+'/masks/'
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+GroupID+'/'


# =============================================================================
# Fichero de parametros
# =============================================================================
# Copiamos el fichero de parametros para asegurarnos de que es el mismo

cmd = 'cp '+results_path+'params.py .'

os.system(cmd)
import params

# Parametros del grupo
bcgID = params.bcgID
z = params.z



# Parametros del instrumento
zp = params.zp
pixelscale = params.pixelscale

# A y Correccion K
k_corrs = params.k_corrs
extinctions = params.extinctions
# sb dimming: -2.5*np.log10(intens) - 10*np.log10(1+z) -- pq el dimming va como (1+z)**-4 -- Calvi+2014  doi: 10.1088/0004-637X/796/2/102
dimming = 10*np.log10(1+z) 


# Parametros de la galaxia
Xc = params.Xc[member]
Yc = params.Yc[member]
ell = params.ell[member]
theta = np.radians(params.theta[member]) # deg. (Desde el eje x, counterclockwise) esto sería nuestro theta
sma0 = params.sma0[member]
maxsma = params.maxsma[member]
step = params.step[member]
mask_file = params.mask_file[member]
os.system('rm params.py')


# =============================================================================
# Leemos la mascara, los datos y los enmascaramos 
# =============================================================================
mask = fits.getdata(masks_path+mask_file)*0
data = fits.getdata(results_path+'PSFcorr-'+band+'.fits')

data_masked  = ma.masked_array(data, mask)






# =============================================================================
# Print check
# =============================================================================
print('# =================================================')
print('               Member: ' +str(member)+'               ')
print('# =================================================')
print('               Fix center: ' +str(fix_center)+'               ')
print('# =================================================')


# =============================================================================
# SB corrections: absorption of our galaxy y k-correction (r band)
# =============================================================================
# Absorcion galactica:  galactic extinction NED 2011ApJ...737..103S  -- 17.7 mag_i
A = extinctions[band] # utilizamos la absorciond de la banda r

k_corr = k_corrs[band]

# =============================================================================
# Elliptical Fitting
# =============================================================================

# Instanciamos nuestra clase
elipses = EllipticalProfile(data_masked, Xc, Yc, sma0, ell, theta,fix_center=fix_center)

# Si no está fijo el centro lo ajustamos                            
# Find the center of a galaxy. If the algorithm is successful the (x, y) coordinates in this EllipseGeometry (i.e. the x0 and y0 attributes) instance will be modified.

if fix_center == False:
    elipses.find_center()


# We perform the elliptical isophote fit with the fit_image method
elipses.fit_image(sma0=sma0, step=step, integrmode='median', sclip=3, nclip=3, fflag=0.1, maxsma=maxsma)

# Ajustamos las zonas externas si se trata de la BCG_sky

#if member=='BCG-sky':
#    elipses.outer_regions(sma_cut=240,step=2*step)

#==============================================================================
# Surface brightness values
#==============================================================================

sbp = elipses.surfbrigth(zp,pixelscale, rms_sky=0, A=A, dimming=dimming, k_corr=k_corr)

# =============================================================================
# Creamos un directorio para guardar los distintos productos del programa
# =============================================================================
try:
    # Create target Directory
    os.mkdir(results_path+member)
except:
    pass

# =============================================================================
# Guardamos la clase instanciada
# =============================================================================
import pickle


pickle.dump(elipses, open(results_path+member+'/'+"main_model.p", "wb" ) )



#==============================================================================
# leemos los datos ajustados
#==============================================================================
r = sbp['sma']
sb = sbp['sb']
err_up = sbp['err up']
err_down = sbp['err down']

ct = sbp['ct']
ct_err = sbp['ct err']


ell = sbp['ell']
ell_err = sbp['ell err']

pa = sbp['pa']
pa_err = sbp['pa err']

xc = sbp['xc']
xc_err = sbp['xc err']
x_median = np.median(xc)

yc = sbp['yc']
yc_err = sbp['yc err']
y_median = np.median(yc)


ell = sbp['ell']
ell_err = sbp['ell err']        
        
pa = sbp['pa']
pa_err = sbp['pa err']

b4 = sbp['b4']
b4_err = sbp['b4 err']

snr = sbp['SNR']

print('X median: ', x_median)
print('Y median: ', y_median)




# =============================================================================
# Figuras de comprobacion        
# =============================================================================
plt.close('all')
#==============================================================================
# Comprobacion: plot de las aperturas
#==============================================================================
plt.figure(1,figsize=(8,8))
plt.title('Group ID: '+GroupID+' - CATAID: '+member,size=25)

# =============================================================================
# Datos sin enmascarar
# =============================================================================
plt.imshow(data, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=10))


# =============================================================================
# Sombreamos las zonas enmascaradas
# =============================================================================
from matplotlib import colors
cmap = colors.ListedColormap(['k'])

shadow = np.nan_to_num(mask,nan=1)
shadow[shadow==0]=np.nan

plt.imshow(shadow,cmap=cmap,  interpolation='nearest',origin='lower',alpha=0.75)

# =============================================================================
# Aperturas
# =============================================================================
elipses.plot(color = 'r')


# =============================================================================
# Ejes en Kpc
# =============================================================================

t_l = baseround(Pix2Kpc(maxsma,scale=pixelscale,z=z))
delta_t = int(t_l/5)

ticks_labels = np.arange(-t_l,t_l+delta_t,delta_t)
ticks = Kpc2Pix(ticks_labels, scale=pixelscale,z=z)


plt.xticks(ticks+Xc,ticks_labels)
plt.yticks(ticks+Yc,ticks_labels)

# =============================================================================
# Limites
# =============================================================================
plt.xlim(Xc-maxsma,Xc+maxsma) 
plt.ylim(Yc-maxsma,Yc+maxsma) 

# =============================================================================
# Etiquetas de los ejes
# =============================================================================

plt.xlabel("kpc",size=20)
plt.ylabel("kpc",size=20)

plt.tight_layout()
plt.savefig(results_path+member+'/MainFitAp.pdf')



# =============================================================================
# Comprobacion: plot de los parametros
# =============================================================================

plt.figure(2, figsize=(8, 8))
plt.suptitle('Group ID: '+GroupID+' - CATAID: '+member,size=25)

plt.subplots_adjust(hspace=0.35, wspace=0.35)

# Ticks y etiquetas para el resto de figuras han de ser los mismo que en la figura 2D pero en la parte positiva

ticks = ticks[[ticks_labels>=0]]
ticks_labels = ticks_labels[ticks_labels>=0]



# =============================================================================
# Ellipticity
# =============================================================================
ax0 = plt.subplot(2, 2, 1) 
plt.plot(r, ell,linestyle='solid')
plt.fill_between(r, ell+ell_err, ell-ell_err, alpha=0.4)
plt.xlim(0,maxsma)
#plt.xscale('symlog',subsx = range(2,10),linthresh=1e-1)
#plt.xticks([0,1,10,100],[0,1,10,100])
plt.xlabel('SMA (pix)',size=14)
plt.ylabel('Ellipticity',size=14)

# =============================================================================
#  Ejes en Kpc
# =============================================================================
ax_twin=ax0.twiny()
ax_twin.tick_params(axis="both",direction="in")
ax_twin.tick_params(axis='x')
ax_twin.set_xlim(0,maxsma)
ax_twin.set_xticks(ticks)
ax_twin.set_xticklabels(ticks_labels)
ax_twin.set_xlabel("SMA (kpc)",size=14)
plt.xlim(0,maxsma)
# =============================================================================


# =============================================================================
# P. Angle
# =============================================================================
ax0 = plt.subplot(2, 2, 2)
plt.plot(r, pa/np.pi*180.,linestyle='solid')
plt.fill_between(r, pa/np.pi*180.+pa_err/np.pi*180., pa/np.pi*180.-pa_err/np.pi*180., alpha=0.4)

plt.xlim(0,maxsma)
#plt.xscale('symlog',subsx = range(2,10))
#plt.xticks([0,1,10,100],[0,1,10,100])
plt.xlabel('SMA (pix)',size=14)
plt.ylabel('PA (deg)',size=14)

# =============================================================================
#  Ejes en Kpc
# =============================================================================
ax_twin=ax0.twiny()
ax_twin.tick_params(axis="both",direction="in")
ax_twin.tick_params(axis='x')
ax_twin.set_xlim(0,maxsma)
ax_twin.set_xticks(ticks)
ax_twin.set_xticklabels(ticks_labels)
ax_twin.set_xlabel("SMA (kpc)",size=14)
plt.xlim(0,maxsma)
# =============================================================================

# =============================================================================
# Xc, Yc 
# =============================================================================
ax0 = plt.subplot(2, 2, 3)
plt.errorbar(r, xc, color='C0', fmt='o',markersize=4,zorder=0)
plt.errorbar(r, yc, color='C0', fmt='o',markersize=4,zorder=0)

plt.axhline(x_median,ls='dashed',color='C1',label='median X0: '+str(round(x_median, 2)),zorder=1)
plt.axhline(y_median,ls='dashed',color='C2',label='median Y0: '+str(round(y_median, 2)),zorder=1)

plt.xlim(0,maxsma)
#plt.xscale('symlog',subsx = range(2,10))
#plt.xticks([0,1,10,100],[0,1,10,100])
plt.xlabel('SMA (pix)',size=14)
plt.ylabel('x0 & Y0',size=14)
plt.legend()

# =============================================================================
#  Ejes en Kpc
# =============================================================================
ax_twin=ax0.twiny()
ax_twin.tick_params(axis="both",direction="in")
ax_twin.tick_params(axis='x')
ax_twin.set_xlim(0,maxsma)
ax_twin.set_xticks(ticks)
ax_twin.set_xticklabels(ticks_labels)
ax_twin.set_xlabel("SMA (kpc)",size=14)
plt.xlim(0,maxsma)
# =============================================================================
# =============================================================================
# B4
# =============================================================================
ax0 = plt.subplot(2, 2, 4)
plt.plot(r, b4, linestyle='solid')
plt.fill_between(r, b4+b4_err, b4-b4_err, alpha=0.4)

plt.axhline(0,ls='dashed',color='grey', zorder=1)

plt.xlim(0,maxsma)
plt.ylim(-0.1, 0.1)

#plt.xscale('symlog',subsx = range(2,10))
#plt.xticks([0,1,10,100],[0,1,10,100])
plt.xlabel('SMA (pix)',size=14)
plt.ylabel('B4',size=14)


# =============================================================================
#  Ejes en Kpc
# =============================================================================
ax_twin=ax0.twiny()
ax_twin.tick_params(axis="both",direction="in")
ax_twin.tick_params(axis='x')
ax_twin.set_xlim(0,maxsma)
ax_twin.set_xticks(ticks)
ax_twin.set_xticklabels(ticks_labels)
ax_twin.set_xlabel("SMA (kpc)",size=14)
plt.xlim(0,maxsma)
# =============================================================================

plt.tight_layout()
plt.savefig(results_path+member+'/MainFitParams.pdf', format='pdf')





#==============================================================================
# PLOT Perfil 
#==============================================================================
plt.close('all')
fig = plt.figure(3,figsize=(10,5))
plt.suptitle('Group ID: '+GroupID+' - CATAID: '+member,size=25)

ax = fig.add_subplot(111)
ax.set_axis_off()
# Perfil en SB
ax1 = fig.add_subplot(121)
ax1.yaxis.set_ticks_position('both') 
ax1.xaxis.set_ticks_position('both') 
ax1.tick_params(direction='in',which='both')  
plt.plot(r, sb,linestyle='solid', label='Mean')
plt.fill_between(r, sb+err_up,sb-err_down, color='lightblue')


# Maximo valor del error que nos creemos
err_up[err_up == 0.0]=err_up.max()
plt.fill_between(r, sb+err_up,sb-err_down, color='k',zorder=0, alpha=0.1,lw=0)

plt.xlim(0, maxsma)
plt.xlabel('SMA (pix)',size=14)

plt.gca().invert_yaxis()
plt.ylabel('$\mu_{'+band+'}$ (mag/arsec$^{2}$)',size=14)



# =============================================================================
#  Ejes en Kpc
# =============================================================================
ax_twin=ax1.twiny()
ax_twin.tick_params(axis="both",direction="in")
ax_twin.tick_params(axis='x')
ax_twin.set_xlim(0,maxsma)
ax_twin.set_xticks(ticks)
ax_twin.set_xticklabels(ticks_labels)
ax_twin.set_xlabel("SMA (kpc)",size=14)
plt.xlim(0,maxsma)
# =============================================================================

# Perfil en cuentas
ax2 = fig.add_subplot(122)
ax2.axhline(0,ls='dashed',color='r')

plt.errorbar(r, ct, yerr=ct_err, fmt='o',markersize=4,zorder=0)
plt.xlabel('SMA (pix)',size=14)
plt.ylabel('Counts',size=14)


# =============================================================================
#  Ejes en Kpc
# =============================================================================
ax_twin=ax2.twiny()
ax_twin.tick_params(axis="both",direction="in")
ax_twin.tick_params(axis='x')
ax_twin.set_xlim(0,maxsma)
ax_twin.set_xticks(ticks)
ax_twin.set_xticklabels(ticks_labels)
ax_twin.set_xlabel("SMA (kpc)",size=14)
plt.xlim(0,maxsma)
# =============================================================================

plt.tight_layout()
plt.show()
plt.savefig(results_path+member+'/MainProfile.pdf')



# =============================================================================
#  Signal to Noise of each ellipse S/N
# =============================================================================
plt.figure(4)
plt.scatter(r, snr)
plt.axhline(1,ls='dashed',color='C1',zorder=1)
plt.xlabel('SMA (pix)')
plt.ylabel('S/N')
plt.ylim(-0.5,10)
plt.tight_layout()
plt.savefig(results_path+member+'/SNR_'+band+'-band.pdf', format='pdf')
# =============================================================================

















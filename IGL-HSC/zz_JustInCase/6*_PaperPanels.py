import numpy as np
import numpy.ma as ma
import os
import pickle
from astropy.table import Table


import matplotlib.pyplot as plt
from matplotlib import gridspec

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

from CML import LogNorm
from CML import Pix2Kpc
from CML import Kpc2Pix

# Division entera en "base"
def baseround(x, base=10):
    return base * round(x/base)

#==============================================================================
# Parametros iniciales
#==============================================================================
GroupID = '400282'
member = 'BCG'


# =============================================================================
# Diccionario de las banda
# =============================================================================
bands = ['g', 'r', 'i']

color = {'g':'C0', 'r':'C2', 'i':'C3'}


# =============================================================================
# SNR limit
# =============================================================================

SNRlim = 0

# =============================================================================
# Fichero de parametros
# =============================================================================
# Copiamos el fichero de parametros para asegurarnos de que es el mismo

cmd = 'cp '+GroupID+'/params.py .'

os.system(cmd)
import params

# Parametros del grupo
bcgID = params.bcgID
z = params.z
dimming = params.dimming


# Parametros del instrumento
zp = params.zp
pixelscale = params.pixelscale

# A , dimming y Correccion K
k_corrs = params.k_corrs
extinctions = params.extinctions

# Parametros de la galaxia
Xc = params.Xc[member]
Yc = params.Yc[member]
step = params.step[member]
mask_file = params.mask_file[member]
os.system('rm params.py')

A = extinctions['r'] # utilizamos la absorciond de la banda r
k_corr = k_corrs['r']




# =============================================================================
# Print check
# =============================================================================
print('# =================================================')
print('               Member: ' +str(member)+'               ')
print('# =================================================')


# =============================================================================
# Leemos la mascara, los datos y los enmascaramo 
# =============================================================================

mask = fits.getdata(GroupID+'/0_masks/'+mask_file)
data = fits.getdata(GroupID+'/0_data/gri-cut.fits')

data_masked  = ma.masked_array(data, mask)

# =============================================================================
# Cargamos el objeto que continene el ajuste original
# =============================================================================
elipses = pickle.load(open(GroupID+'/'+member+'/'+"main_model.p", "rb" ) )



# =============================================================================
# Calculamos el maximo sma considerado en entre las 3 bandas
# =============================================================================
maxsma = []

for band in bands:
    path = GroupID+'/'+member+'/sbp_'+band+'-band.fits'
    sbp = Table.read(path)
    
    index = sbp['SNR'] >= SNRlim
    
    sb = sbp['sb'][index]
    maxsma.append(np.max(sbp['sma'][index]))

maxsma = np.max(maxsma)



# =============================================================================
# Ticks
# =============================================================================
# Kpc

t_l = baseround(Pix2Kpc(maxsma,scale=pixelscale,z=z))
delta_t = int(t_l/5)

kpc_ticks_labels = np.arange(0,t_l,delta_t)
kpc_ticks = Kpc2Pix(kpc_ticks_labels, scale=pixelscale,z=z)

# arcsec

t_l = baseround(maxsma*pixelscale)
delta_t = int(t_l/5)

arcsec_ticks_labels =  np.arange(0,t_l,delta_t)
arcsec_ticks = arcsec_ticks_labels/pixelscale


# =============================================================================
# Comenzamos la figura        
# =============================================================================
plt.close('all')
fig = plt.figure(1,figsize=(12.5,6.5))

# Grid de subplots
gs = gridspec.GridSpec(3, 6)



#==============================================================================
# PLOT Perfil 
#==============================================================================

# =============================================================================
band = 'g'    
# =============================================================================

ColR = fig.add_subplot(gs[0:1,-2:])

path = GroupID+'/'+member+'/sbp_'+band+'-band.fits'
sbp = Table.read(path)


# =============================================================================
#  Cargamos los datos y consideramos solos los puntos con SNR mayor que SNRlim
# =============================================================================
index = sbp['SNR'] >= SNRlim

sb = sbp['sb'][index]
r = sbp['sma'][index]
err_up = sbp['err up'][index]
err_down = sbp['err down'][index]


ColR.plot(r, sb,linestyle='solid', label='Mean', color = color[band])
plt.fill_between(r, sb-err_down,sb+err_up, alpha=0.4,color = color[band],lw=0)


ColR.yaxis.set_label_position("right")
ColR.yaxis.tick_right()
ColR.tick_params(axis="both",direction="in")


plt.gca().invert_yaxis()

plt.ylabel('$\mu_{'+band+'}$ (mag/arsec$^{2}$)',size=14)


# =============================================================================
#  Ejes en Kpc
# =============================================================================

plt.xticks(kpc_ticks,[])
plt.grid(ls=':')
plt.xlim(0,maxsma)
# =============================================================================
#  Ejes en arcsec
# =============================================================================

ax_twin=ColR.twiny()
ax_twin.tick_params(axis="both",direction="in")
ax_twin.set_xticks(arcsec_ticks)
ax_twin.set_xticklabels(arcsec_ticks_labels)
ax_twin.set_xlabel("SMA (arcsec)",size=12)
ax_twin.set_xlim(0,maxsma)



# =============================================================================
band = 'r'    
# =============================================================================

ColR = fig.add_subplot(gs[1:2,-2:])

path = GroupID+'/'+member+'/sbp_'+band+'-band.fits'
sbp = Table.read(path)


# =============================================================================
#  Cargamos los datos y consideramos solos los puntos con SNR mayor que SNRlim
# =============================================================================
index = sbp['SNR'] >= SNRlim

sb = sbp['sb'][index]
r = sbp['sma'][index]
err_up = sbp['err up'][index]
err_down = sbp['err down'][index]


ColR.plot(r, sb,linestyle='solid', label='Mean', color = color[band])
plt.fill_between(r, sb-err_down,sb+err_up, alpha=0.4,color = color[band],lw=0)




ColR.yaxis.set_label_position("right")
ColR.yaxis.tick_right()
ColR.tick_params(axis="both",direction="in")


plt.gca().invert_yaxis()

plt.ylabel('$\mu_{'+band+'}$ (mag/arsec$^{2}$)',size=14)


# =============================================================================
#  Ejes en Kpc
# =============================================================================

plt.xticks(kpc_ticks,[])
plt.grid(ls=':')
plt.xlim(0,maxsma)
# =============================================================================
#  Ejes en arcsec
# =============================================================================

ax_twin=ColR.twiny()
ax_twin.tick_params(axis="both",direction="in")
ax_twin.set_xticks(arcsec_ticks)
ax_twin.set_xticklabels([])
ax_twin.set_xlim(0,maxsma)


# =============================================================================
band = 'i'    
# =============================================================================

ColR = fig.add_subplot(gs[2:,-2:])

path = GroupID+'/'+member+'/sbp_'+band+'-band.fits'
sbp = Table.read(path)


# =============================================================================
#  Cargamos los datos y consideramos solos los puntos con SNR mayor que SNRlim
# =============================================================================
index = sbp['SNR'] >= SNRlim

sb = sbp['sb'][index]
r = sbp['sma'][index]
err_up = sbp['err up'][index]
err_down = sbp['err down'][index]


ColR.plot(r, sb,linestyle='solid', label='Mean', color = color[band])
plt.fill_between(r, sb-err_down,sb+err_up, alpha=0.4,color = color[band],lw=0)

ColR.yaxis.set_label_position("right")
ColR.yaxis.tick_right()
ColR.tick_params(axis="both",direction="in")


plt.gca().invert_yaxis()

plt.ylabel('$\mu_{'+band+'}$ (mag/arsec$^{2}$)',size=14)


# =============================================================================
#  Ejes en Kpc
# =============================================================================

plt.xlabel("SMA (kpc)",size=12)
plt.grid(ls=':')
plt.xticks(kpc_ticks,kpc_ticks_labels)
plt.xlim(0,maxsma)
# =============================================================================
#  Ejes en arcsec
# =============================================================================

ax_twin=ColR.twiny()
ax_twin.tick_params(axis="both",direction="in")
ax_twin.set_xticks(arcsec_ticks)
ax_twin.set_xticklabels([])
ax_twin.set_xlim(0,maxsma)



# =============================================================================
# 2-D
# =============================================================================

CentFig = fig.add_subplot(gs[:,1:4])
CentFig.yaxis.set_label_coords(-0.04,0.5)
CentFig.tick_params(axis="both",direction="inout",width=1,length=7)


plt.title('Group ID: '+GroupID+' - Member ID: '+member,size=17,y=1.08)


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

CentFig.imshow(shadow,cmap=cmap,  interpolation='nearest',origin='lower',alpha=0.75,aspect='auto')

# =============================================================================
# Aperturas
# =============================================================================
elipses.plot(color = 'w',smamax=maxsma)


# =============================================================================
# Ejes en Kpc
# =============================================================================

kpc_ticks_2d = np.concatenate([-kpc_ticks[1:],kpc_ticks])
kpc_ticks_labels_2d = np.concatenate([-kpc_ticks_labels[1:],kpc_ticks_labels])

CentFig.set_xticks(kpc_ticks_2d+Xc)
CentFig.set_xticklabels(kpc_ticks_labels_2d,size=10)
CentFig.set_yticks(kpc_ticks_2d+Yc)
CentFig.set_yticklabels(kpc_ticks_labels_2d,size=10)

CentFig.set_xlabel("kpc",size=15)
CentFig.set_ylabel("kpc",size=15)

CentFig.set_xlim(Xc-maxsma,Xc+maxsma) 
CentFig.set_ylim(Yc-maxsma,Yc+maxsma) 

# =============================================================================
#  Ejes en arcsec
# =============================================================================

arcsec_ticks_2d = np.concatenate([-arcsec_ticks[1:],arcsec_ticks])
arcsec_ticks_labels_2d = np.concatenate([-arcsec_ticks_labels[1:],arcsec_ticks_labels])

# Eje X

ax_twin=CentFig.twiny()
ax_twin.tick_params(axis="both",direction="in")

ax_twin.set_xticks(arcsec_ticks_2d+Xc)
ax_twin.set_xticklabels(arcsec_ticks_labels_2d)

ax_twin.set_xlim(Xc-maxsma,Xc+maxsma) 
ax_twin.set_xlabel("arcsec",size=15) 

# Eje Y

ay_twin=CentFig.twinx()
ay_twin.tick_params(axis="both",direction="in")

ay_twin.set_yticks(arcsec_ticks_2d+Xc)
ay_twin.set_yticklabels(arcsec_ticks_labels_2d)

ay_twin.set_ylim(Xc-maxsma,Xc+maxsma) 
ay_twin.set_ylabel("arcsec",size=15,rotation=270) 
ay_twin.yaxis.set_label_coords(1.06,0.5)




# =============================================================================
# Fit Parameters
# =============================================================================


#==============================================================================
# Surface brightness values
#==============================================================================

sbp = elipses.surfbrigth(zp, pixscale=pixelscale, rms_sky=0, A=A, dimming=dimming, k_corr=k_corr)


#==============================================================================
# leemos los datos ajustados
#==============================================================================
r = sbp['sma']

ell = sbp['ell']
ell_err = sbp['ell err']


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



# =============================================================================
# Ellipticity
# =============================================================================

ColIz0 = fig.add_subplot(gs[0:1,0:1])
ColIz0.tick_params(axis="both",direction="in")



plt.plot(r, ell,linestyle='solid',color='C1')
plt.fill_between(r, ell+ell_err, ell-ell_err, alpha=0.4,color='C1')

# =============================================================================
#  Ejes en Kpc
# =============================================================================

plt.xticks(kpc_ticks,[])
plt.grid(ls=':')
plt.xlim(0,maxsma)

plt.ylabel('Ellipticity',size=14)


# =============================================================================
#  Ejes en arcsec
# =============================================================================

ax_twin=ColIz0.twiny()
ax_twin.tick_params(axis="both",direction="in")
ax_twin.set_xticks(arcsec_ticks)
ax_twin.set_xticklabels(arcsec_ticks_labels)
ax_twin.set_xlabel("SMA (arcsec)",size=12)
ax_twin.set_xlim(0,maxsma)



# =============================================================================
# P. Angle
# =============================================================================
ColIz1 = fig.add_subplot(gs[1:2,0:1])
ColIz1.tick_params(axis="both",direction="in")


plt.plot(r, pa/np.pi*180.,linestyle='solid',color='C1')
plt.fill_between(r, pa/np.pi*180.+pa_err/np.pi*180., pa/np.pi*180.-pa_err/np.pi*180., alpha=0.4,color='C1')

# =============================================================================
#  Ejes en Kpc
# =============================================================================

plt.xticks(kpc_ticks,[])
plt.grid(ls=':')
plt.xlim(0,maxsma)

plt.ylabel('PA',size=14)


# =============================================================================
#  Ejes en arcsec
# =============================================================================

ax_twin=ColIz1.twiny()
ax_twin.tick_params(axis="both",direction="in")
ax_twin.set_xticks(arcsec_ticks)
ax_twin.set_xticklabels([])
ax_twin.set_xlim(0,maxsma)





# =============================================================================
# =============================================================================
# B4
# =============================================================================
ColIz2 = fig.add_subplot(gs[2:,0:1])
ColIz2.tick_params(axis="both",direction="in")

ColIz2.yaxis.set_label_coords(-0.2,0.5)



plt.plot(r, b4, linestyle='solid',color='C1')
plt.fill_between(r, b4+b4_err, b4-b4_err, alpha=0.4,color='C1')

plt.axhline(0,ls='dashed',color='grey', zorder=1)

# =============================================================================
#  Ejes en Kpc
# =============================================================================

plt.xticks(kpc_ticks,kpc_ticks_labels)
plt.grid(ls=':')
plt.xlim(0,maxsma)

plt.xlabel('SMA (kpc)',size=14)
plt.ylabel('B4',size=14)


# =============================================================================
#  Ejes en arcsec
# =============================================================================

ax_twin=ColIz2.twiny()
ax_twin.tick_params(axis="both",direction="in")
ax_twin.set_xticks(arcsec_ticks)
ax_twin.set_xticklabels([])
ax_twin.set_xlim(0,maxsma)


plt.subplots_adjust(wspace=0.4, hspace=0.1,bottom=0.085,left=0.055,right=0.95,top=0.88)    				

plt.savefig(GroupID+'/'+member+'/PaperFig-'+member+'.pdf')










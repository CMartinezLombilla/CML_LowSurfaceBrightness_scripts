import numpy as np
import numpy.ma as ma
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.close('all')
plt.ioff()

plt.rcParams["font.family"] = "serif"
plt.rc('text', usetex=False)

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord


# =============================================================================
# importamos CML
# =============================================================================

import sys
sys.path.append('/Users/felipe/Dropbox/CML/') # path del CML Feli
sys.path.append('/Users/z3530031/Dropbox/CML/') # path del CML Cris

sys.path.append('/Users/felipe/Dropbox/CML/WHT/') # path del fichero de parametros Feli
sys.path.append('/Users/z3530031/Dropbox/CML/WHT/') # path del fichero de parametros Cris

from CML import LogNorm
from CML import CircularApProfile

# Division entera en "base"
def baseround(x, base=10):
    return base * round(x/base)


# =============================================================================
# Parametros iniciales del grupo/imagen/IGL/aperturas
# =============================================================================
groupID = '400138'
member = 'IGL'
band = 'i'
sigma = 3
nbins = 20


# =============================================================================
# Paths
# =============================================================================
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'


# =============================================================================
# Creamos un directorio para guardar los distintos productos del programa
# =============================================================================
try:
    # Create target Directory
    os.mkdir(results_path+'IGL')
except:
    pass


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
pixscale = params.pixelscale

# A y Correccion K
k_corrs = params.k_corrs
extinctions = params.extinctions
SkyRMS = params.SkyRMS

# sb dimming: -2.5*np.log10(intens) - 10*np.log10(1+z) -- pq el dimming va como (1+z)**-4 -- Calvi+2014  doi: 10.1088/0004-637X/796/2/102
dimming = 10*np.log10(1+z) 

# Parametros de la galaxia
Xc = params.Xc[member]
Yc = params.Yc[member]
sma0 = params.sma0[member]
maxsma = params.maxsma[member]
step = params.step[member]
mask_file = params.mask_file[member]
os.system('rm params.py')


# =============================================================================
# Leemos la mascara, los datos y los enmascaramo 
# =============================================================================
data = fits.getdata(results_path+'IGLlight-'+band+'.fits')
mask = fits.getdata(results_path+'masks/'+mask_file)
data_masked  = ma.masked_array(data, mask)



# =============================================================================
# Aperturas
# =============================================================================
ap = CircularApProfile(data_masked, Xc, Yc, nbins=nbins, r_min=sma0, npix=maxsma)


# =============================================================================
# Perfil
# =============================================================================

t = ap.sbprofile(zp, pixscale, A=extinctions[band],rms_sky= SkyRMS[band], dimming=dimming, k_corr=k_corrs[band], sigma=3)

r = t['r'][:-2]
sb = t['sb'][:-2]
err_up = t['err up'][:-2]
err_down = t['err down'][:-2]





# =============================================================================
# Figuras de comprobacion        
# =============================================================================
plt.close('all')
# =============================================================================
# PDF para las imagenes de control
# =============================================================================
pdf = PdfPages(results_path+'IGL/IGL_profile.pdf')

               
#==============================================================================
# Comprobacion: plot de las aperturas
#==============================================================================
plt.figure(1,figsize=(8,8))
plt.title('Group ID: '+groupID,size=25)
# Datos sin enmascarar
plt.imshow(data, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=10))
ap.plot(color='white', ls=':')
# =============================================================================
# Sombreamos las zonas enmascaradas
# =============================================================================
from matplotlib import colors
cmap = colors.ListedColormap(['k'])

shadow = np.nan_to_num(mask,nan=1)
shadow[shadow==0]=np.nan

plt.imshow(shadow,cmap=cmap,  interpolation='nearest',origin='lower',alpha=0.7)
plt.xlim(np.shape(data_masked)[0]/2-maxsma*1.2, np.shape(data_masked)[0]/2+maxsma*1.2)
plt.ylim(np.shape(data_masked)[0]/2-maxsma*1.2, np.shape(data_masked)[0]/2+maxsma*1.2)



# =============================================================================
# Perfil sb
# =============================================================================
max_x = 240
# =============================================================================
# Ticks en Kpc
# =============================================================================
from CML import Pix2Kpc
from CML import Kpc2Pix
t_l = baseround(Pix2Kpc(max_x,scale=pixscale,z=z))
delta_t = int(t_l/5)

ticks_labels = np.arange(0,t_l+delta_t,delta_t)
ticks = Kpc2Pix(ticks_labels, scale=pixscale,z=z)



fig,ax = plt.subplots(figsize=(8, 5))


ax.errorbar(r, sb, yerr=[err_down,err_up], marker='.',mec='k',mew=0.2, ms=7,linestyle='-',lw=1.3,color='r')
ax.grid('dashed')
ax.set_ylim(31.5,23.9)
ax.set_xlim(0,max_x)

ax.yaxis.set_ticks_position('both') 
ax.xaxis.set_ticks_position('both') 
ax.set_xlabel('SMA (pix)',size=14)

ax.set_ylabel('$\mu_{'+band+'}$ (mag/arsec$^{2}$)',size=14)
    
    
# =============================================================================
#  Ejes en Kpc
# =============================================================================
ax_twin=ax.twiny()
ax_twin.tick_params(axis="both",direction="in")
ax_twin.tick_params(axis='x')
ax_twin.set_xlim(0,max_x)
ax_twin.set_xticks(ticks)
ax_twin.set_xticklabels(ticks_labels)
ax_twin.set_xlabel("SMA (kpc)",size=14)
plt.xlim(0,max_x)


plt.show()


for fig in range(1, plt.figure().number): ## will open an empty extra figure :(
    pdf.savefig( fig )
pdf.close()














from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import numpy as np
import os
from astropy.table import Table

plt.rcParams["font.family"] = "serif"

# =============================================================================
# importamos CML
# =============================================================================

import sys
sys.path.append('/Users/felipe/Dropbox/CML/') # path del CML Feli
sys.path.append('/Users/z3530031/Dropbox/CML/') # path del CML Cris

sys.path.append('/Users/felipe/Dropbox/CML/WHT/') # path del fichero de parametros Feli
sys.path.append('/Users/z3530031/Dropbox/CML/WHT/') # path del fichero de parametros Cris

from CML import Pix2Kpc
from CML import Kpc2Pix
# Division entera en "base"
def baseround(x, base=10):
    return base * round(x/base)



#==============================================================================
#Initial parameters
#==============================================================================

GroupID = '400282'
member = 'BCG'
SNR_cut = -10
# =============================================================================
# Diccionarios para las bandas
# =============================================================================

bands = ['g', 'r', 'i']

color = {'g':'C0', 'r':'C2', 'i':'C3'}

yticks = {'g':np.int_(np.arange(18,34,2)),
          'r':np.int_(np.arange(16,32,2)),
          'i':np.int_(np.arange(16,32,2))}

ylims = {'g':(17.8,33),
          'r':(15.85,31.8),
          'i':(15.85,31.8)}


# Diccionario para los subplots 
plt.close('all')
fig = plt.figure(1,figsize=(8,8))





# =============================================================================
# Fichero de parametros 
# =============================================================================

cmd = 'cp '+GroupID+'/params.py .'

os.system(cmd)
import params

# Parametros del grupo
bcgID = params.bcgID
z = params.z
zp = params.zp
maxsma = params.maxsma[member]
pixelscale = params.pixelscale

sb_limit = params.sb_limit

os.system('rm params.py')

# =============================================================================
# Print check
# =============================================================================
print('# =================================================')
print('               Member: ' +str(member)+'               ')
print('# =================================================')


# =============================================================================
# Ticks en Kpc
# =============================================================================

t_l = baseround(Pix2Kpc(maxsma,scale=pixelscale,z=z))
delta_t = int(t_l/5)

ticks_labels = np.arange(0,t_l+delta_t,delta_t)
ticks = Kpc2Pix(ticks_labels, scale=pixelscale,z=z)




# =============================================================================
# Titulos en el subplot vacio
# =============================================================================
ax = fig.add_subplot(221)

ax.axis('off')
ax.annotate('Goup ID: '+GroupID+'\n\n  Member: '+member,(0.5, 0.5), xycoords='axes fraction', ha='center', va='center', size = 22)

# =============================================================================
# Loop sobre las bandas
# =============================================================================

for band in bands:
    
    path = GroupID+'/'+member+'/sbp_'+band+'-band.fits'
    
    sbp = Table.read(path)

    # =============================================================================
    # Leemos todas las variables
    # =============================================================================
    # Nos quedamos con el los puntos con SNR mayor que 1
    
    index = sbp['SNR']>= SNR_cut
    
    sb_uncorrected = - 2.5*np.log10(sbp['ct'][index]) + zp + 5*np.log10(pixelscale)


    sb = sbp['sb'][index]
    r = sbp['sma'][index]
    err_up = sbp['err up'][index]
    err_down = sbp['err down'][index]
    
    maxsma = r.max()
    
    # =============================================================================
    # PLot de los perfiles
    # =============================================================================
    ind = bands.index(band)
    ax = fig.add_subplot(222+ind)
    plt.title(band+'-band', y=0.85,x=0.85,size=15)
    ax.yaxis.set_ticks_position('both') 
    ax.xaxis.set_ticks_position('both') 
    ax.tick_params(direction='in',which='both')  
        
    #ax.errorbar(r, sb, yerr=[err_down,err_up], marker='.',mec='k',mew=0.2, ms=7,linestyle='-',lw=1.3,color=color[band])

    ax.plot(r, sb,linestyle='solid', color = color[band])
    plt.fill_between(r, sb-err_down,sb+err_up, alpha=0.2,color = color[band],lw=0)
    
    # SBP uncorrected
    #ax.plot(r, sb_uncorrected, ls='dashed',lw=0.8, color = 'grey')
    
    # Limiting SB 
    #ax.axhline(sb_limit[band],color='k',lw=0.9,ls=':',label='$\mu_{lim}$: '+str(sb_limit[band]))
    #ax.legend(loc=(0.05,0.85),frameon=False)
    
    # Maximo valor del error que nos creemos
    err_up[err_up == 0.0]=err_up.max()
    plt.fill_between(r, sb+err_up,sb-err_down, color='k',zorder=0, alpha=0.1,lw=0)

    
    
    plt.xlim(0, maxsma)
    #plt.xscale('symlog',subs = range(2,10))
    plt.xlabel('SMA (pix)',size=14)
    
    plt.gca().invert_yaxis()
    plt.ylabel('$\mu_{'+band+'}$ (mag/arsec$^{2}$)',size=14)
    
    
    # =============================================================================
    #  Ejes en Kpc
    # =============================================================================
    ax_twin=ax.twiny()
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

plt.savefig(GroupID+'/'+member+'/AllBandProfiles.pdf')

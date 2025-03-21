import numpy as np
import numpy.ma as ma
import os
import pickle

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.close('all')
plt.ioff()

plt.rcParams["font.family"] = "serif"
plt.rc('text', usetex=False)

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.table import Table


# =============================================================================
# importamos CML
# =============================================================================

import sys
sys.path.append('/Users/felipe/Dropbox/CML/') # path del CML Feli
sys.path.append('/Users/z3530031/Dropbox/CML/') # path del CML Cris

sys.path.append('/Users/felipe/Dropbox/CML/WHT/') # path del fichero de parametros Feli
sys.path.append('/Users/z3530031/Dropbox/CML/WHT/') # path del fichero de parametros Cris

from CML import LogNorm
from CML import EllipticalApProfile
from CML import isophotes_prof
from CML import save_fits_image
from CML import Counts2Sb
from CML import com
from CML import dist_galaxies
from CML import Pix2Kpc


# Division entera en "base"
def baseround(x, base=10):
    return base * round(x/base)


# =============================================================================
# Parametros iniciales del grupo/imagen/IGL/aperturas
# =============================================================================
groupID = '400138'
members_core = ['BCG', '1660730', '1660646', '1660615',]
members = ['BCG', '1660646','1660730',]
bands = ['g', 'r', 'i']

# Params profile
sma0 = 2
nbins = 27

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

try:
    # Create target Directory
    os.mkdir(results_path+'Isophotes/IGL/profiles_isoph')
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
SkyRMS = params.SkyRMS

# sb dimming: -2.5*np.log10(intens) - 10*np.log10(1+z) -- pq el dimming va como (1+z)**-4 -- Calvi+2014  doi: 10.1088/0004-637X/796/2/102
dimming = 10*np.log10(1+z) 

mask_file = params.mask_file['IGL']
   
#Posiciones centrales del core
Xcore = []
Ycore = []
for memberc in members_core:
    Xcore.append(int(round(params.Xc[memberc],0)))
    Ycore.append(int(round(params.Yc[memberc],0)))    

# Galaxies params
x_gals = []
y_gals = []
ell_gals = []
pa_gals = [] # PA 2d-imfit + 90 deg


for member in members:              
    # Galaxy params
    x_gals.append(params.Xc[member])
    y_gals.append(params.Yc[member])
    ell_gals.append(1-params.ell[member])
    pa_gals.append(params.theta[member]) # deg. (Desde el eje x, counterclockwise) esto sería nuestro theta

os.system('rm params.py')


# Find the center of 'light' of the core region
Gmasked = fits.getdata(results_path+'masks/GroupMasked-i.fits')
Xc_core, Yc_core = com(Gmasked, Xcore, Ycore)



# =============================================================================
# Obtain sbp using distances map
# =============================================================================


for structure in ['IGL', 'group']:
  
    for band in bands:          
                
        # Dictionaries for group and IGL-only data
        isophotes = fits.getdata(results_path+'Isophotes/'+'isophotes_'+band+'.fits')

        data = {'IGL': results_path+'Dist_IGLlight-'+band+'.fits', 
                'group': results_path+'PSF/PSFcorr-'+band+'-SkySubDist.fits'}
        
        mask = {'IGL': results_path+'masks/'+mask_file+'-'+band+'.fits', 
                'group': results_path+'masks/CoreMask-'+band+'.fits'}
              
        # Masking the data    
        im  = np.copy(fits.getdata(data[structure]))
        im[np.isnan(fits.getdata(mask[structure]))] = np.nan
        
        # =============================================================================
        
        # Get the profiles
        prof = isophotes_prof(im, isophotes, Xc_core, Yc_core, SkyRMS[band], zp, pixscale, A=extinctions[band], dimming=dimming, k_corr=k_corrs[band])
        prof.write(results_path+'Isophotes/IGL/profiles_isoph/'+structure+'_sbp_'+band+'.fits', overwrite=True)
    
        #r = prof['r']
        sb = prof['sb']
        err_up = prof['err up']
        err_down = prof['err down']
        snr = prof['SNR']
       
        
        # Get the distances map               
        r = dist_galaxies(isophotes, x_gals, y_gals, ell_gals, pa_gals)
        #r = Pix2Kpc(r, pixscale, z=z)
        
        '''
        # Get the profiles with No Corrections
        prof_NC = dist_ellipse_ct_prof(im, dist_map, sma0, step, nbins, SkyRMS[band],zp, pixscale, A=0, dimming=0, k_corr=0)
        prof_NC.write(results_path+'IGL/profiles/dist_'+structure+'_sbp_NC_'+band+'.fits', overwrite=True)
    
        r_NC = prof_NC['r']
        sb_NC = prof_NC['sb']
        err_up_NC = prof_NC['err up']
        err_down_NC = prof_NC['err down']
        snr_NC = prof_NC['SNR']
        '''

        # =============================================================================
        # Figures
        # =============================================================================
        plt.close('all')
        pdf = PdfPages(results_path+'Isophotes/IGL/profiles_isoph/'+structure+'_profile-'+band+'.pdf')
               
        
        # =============================================================================
        # Perfil sb
        # =============================================================================
        #max_x = 300
        # =============================================================================
        # Ticks en Kpc
        # =============================================================================
        from CML import Pix2Kpc
        from CML import Kpc2Pix
        
        fig,ax = plt.subplots(figsize=(8, 5))
        
        ax.errorbar(r, sb, yerr=[err_down,err_up], marker='.',mec='k',mew=0.2, ms=7,linestyle='-',lw=1.3,color='r', label='IGL SBP')
        #ax.errorbar(r_NC, sb_NC, yerr=[err_down_NC,err_up_NC], marker='.',mec='k',mew=0.2, ms=7,linestyle=':',lw=1, color='r', alpha=0.6, label='IGL SBP No Corr')
    
        ax.grid('dashed')
        ax.legend()
        ax.set_ylim(32.9,17.5)
        #ax.set_xlim(0,max_x)
        
        ax.yaxis.set_ticks_position('both') 
        ax.xaxis.set_ticks_position('both') 
        ax.set_xlabel('Isophote as in fits image',size=14)
        
        ax.set_ylabel('$\mu_{'+band+'}$ (mag/arsec$^{2}$)',size=14)
            
        '''    
        # =============================================================================
        #  Ejes en Kpc
        # =============================================================================
        t_l = baseround(Pix2Kpc(max_x,scale=pixscale,z=z))
        delta_t = int(t_l/5)
        
        ticks_labels = np.arange(0,t_l+delta_t,delta_t)
        ticks = Kpc2Pix(ticks_labels, scale=pixscale,z=z)
        
        ax_twin=ax.twiny()
        ax_twin.tick_params(axis="both",direction="in")
        ax_twin.tick_params(axis='x')
        ax_twin.set_xlim(0,max_x)
        ax_twin.set_xticks(ticks)
        ax_twin.set_xticklabels(ticks_labels)
        ax_twin.set_xlabel("SMA (kpc)",size=14)
        plt.xlim(0,max_x)
        '''
        
        
        # =============================================================================
        #  Signal to Noise of each ellipse S/N
        # =============================================================================
        plt.figure(2)
        plt.scatter(r, snr)
        plt.axhline(1,ls='dashed',color='C1',zorder=1)
        plt.xlabel('Isophote as in fits image')
        plt.ylabel('S/N')
        plt.ylim(-0.5,200)
        plt.tight_layout()
        # =============================================================================  
    
        
        plt.show()
        
        
    for fig in range(1, plt.figure().number): ## will open an empty extra figure :(
        pdf.savefig( fig )
    pdf.close()
        
        
        










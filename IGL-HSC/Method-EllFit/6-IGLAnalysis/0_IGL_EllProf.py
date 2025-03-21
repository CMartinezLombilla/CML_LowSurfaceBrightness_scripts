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


# Division entera en "base"
def baseround(x, base=10):
    return base * round(x/base)


# =============================================================================
# Parametros iniciales del grupo/imagen/IGL/aperturas
# =============================================================================
groupID = '400138'
members = ['BCG', '1660730', '1660646']
bands = ['g', 'r', 'i']
fitting = 'custom' # 'custom' or 'BCGsky'
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

try:
    # Create target Directory
    os.mkdir(results_path+'IGL/profiles')
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


# Parametros
sma0 = params.sma0['IGL']
maxsma = params.maxsma['IGL']
step = params.step['IGL']
mask_file = params.mask_file['IGL']
ell = params.ell['IGL']
theta = params.theta['IGL']


Xcoords = []
Ycoords = []
for member in members:
    Xcoords.append(int(round(params.Xc[member],0)))
    Ycoords.append(int(round(params.Yc[member],0)))
    
os.system('rm params.py')





if fitting == 'custom':
    # =============================================================================
    # Usamos el ell y pa de las ellipses en BCGsky y diferentes steps/bin sizes 
    # =============================================================================
    ell_sky = pickle.load(open(results_path+'BCGsky/main_model.p', 'rb' ) )
      
    sma_arr= np.logspace(np.log10(sma0),np.log10(maxsma),nbins)
    
    ell_arr = np.zeros_like(sma_arr)+np.median(ell_sky.isolist.eps)
    pa_arr = np.zeros_like(sma_arr)+np.median(ell_sky.isolist.pa)
    
    # =============================================================================

else:
    # =============================================================================
    # Usamos el ajuste original de ellipses en BCGsky
    # =============================================================================
    ell_sky = pickle.load(open(results_path+'BCGsky/main_model.p', 'rb' ) )
    
    sma_arr=ell_sky.isolist.sma[-nbins:]
    ell_arr=ell_sky.isolist.eps[-nbins:] 
    pa_arr=ell_sky.isolist.pa[-nbins:]
    # =============================================================================


for structure in ['IGL', 'group']:
    
    for band in bands:
        
        # =============================================================================
        # Find the center of masss of the system
        # =============================================================================
        Gmasked = fits.getdata(results_path+'masks/GroupMasked-'+band+'.fits')
        weights_list = []
        
        for ii in np.arange(len(members)):
            weights_list.append(Gmasked[Ycoords[ii], Xcoords[ii]])
        
        weights = np.array(weights_list)
        
        com_x = sum(Xcoords * weights) / sum(weights)
        com_y = sum(Ycoords * weights) / sum(weights)
        
        # Center of mass of the system (use it for the extraction of the IGL)
        Xc = com_x 
        Yc = com_y
        
        
        
        # =============================================================================
        # Dictionaries for group and IGL-only data, we masl them too
        # =============================================================================
        
        data = {'IGL': results_path+'IGLlight-'+band+'.fits', 
                'group': results_path+'PSF/PSFcorr-'+band+'-SkySub.fits'}
        
        mask = {'IGL': results_path+'masks/'+mask_file+'-'+band+'.fits', 
                'group': results_path+'masks/CoreMask-'+band+'.fits'}
        
        data = fits.getdata(data[structure])
        mask = fits.getdata(mask[structure])
        data_masked  = ma.masked_array(data, mask)
     
        
        # =============================================================================
        # Elliptical Fitting - usando las aperturas especificadas antes
        # =============================================================================
        # Instanciamos nuestra clase
        elipses = EllipticalApProfile(data_masked, Xc, Yc, sma_arr=sma_arr, ell_arr=ell_arr, pa_arr=pa_arr)
    
        # Plot de las aperturas
        #elipses.plot()
       
        
        #==============================================================================
        # Surface brightness values
        #==============================================================================    
        sbp_NC = elipses.surfbrigth(zp, pixscale, rms_sky=SkyRMS[band], A=0, dimming=0, k_corr=0)
        
        r_NC = sbp_NC['r'][:-1]
        sb_NC = sbp_NC['sb'][:-1]
        err_up_NC = sbp_NC['err up'][:-1]
        err_down_NC = sbp_NC['err down'][:-1]
        snr_NC = sbp_NC['SNR']
        
        sbp = elipses.surfbrigth(zp, pixscale, rms_sky=SkyRMS[band], A=extinctions[band], dimming=dimming, k_corr=k_corrs[band])
        
        sbp.write(results_path+'IGL/profiles/'+structure+'_sbp_'+band+'.fits', overwrite=True)
        
        r = sbp['r'][:-1]
        sb = sbp['sb'][:-1]
        err_up = sbp['err up'][:-1]
        err_down = sbp['err down'][:-1]
        snr = sbp['SNR']
        
        # =============================================================================
        # Guardamos la clase instanciada
        # =============================================================================
        pickle.dump(elipses, open(results_path+'IGL/profiles/main_model-'+structure+'-'+band+'.p', 'wb' ) )
        
        
        
        # =============================================================================
        # Figuras de comprobacion        
        # =============================================================================
        plt.close('all')
        # =============================================================================
        # PDF para las imagenes de control
        # =============================================================================
        pdf = PdfPages(results_path+'IGL/profiles/'+structure+'_profile-'+band+'.pdf')     
                       
        #==============================================================================
        # Comprobacion: plot de las aperturas
        #==============================================================================
        plt.figure(1,figsize=(8,8))
        plt.title('Group ID: '+groupID,size=25)
        # Datos sin enmascarar
        plt.imshow(data, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=10))
        elipses.plot(color='white', ls=':')
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
        max_x = 580
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
        
        
        ax.errorbar(r, sb, yerr=[err_down,err_up], marker='.',mec='k',mew=0.2, ms=7,linestyle='-',lw=1.3,color='r', label='IGL SBP')
        ax.errorbar(r_NC, sb_NC, yerr=[err_down_NC,err_up_NC], marker='.',mec='k',mew=0.2, ms=7,linestyle=':',lw=1, color='r', alpha=0.6, label='IGL SBP No Corr')
    
        ax.grid('dashed')
        ax.legend()
        ax.set_ylim(32.9,22.5)
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
        
        
        
        # =============================================================================
        #  Signal to Noise of each ellipse S/N
        # =============================================================================
        plt.figure(3)
        plt.scatter(r, snr[:-1])
        plt.axhline(1,ls='dashed',color='C1',zorder=1)
        plt.xlabel('SMA (pix)')
        plt.ylabel('S/N')
        plt.ylim(-0.5,200)
        plt.tight_layout()
        # =============================================================================  
    
                
        
        for fig in range(1, plt.figure().number): ## will open an empty extra figure :(
            pdf.savefig( fig )
        pdf.close()
        
    
    











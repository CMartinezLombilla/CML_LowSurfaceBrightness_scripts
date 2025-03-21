import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"

# =============================================================================
# import CML
# =============================================================================

import sys
sys.path.append('/Users/z3530031/Dropbox/CML/') # path del CML Cris

from CML import LogNorm
from CML import save_fits_image


# =============================================================================
# Parametros iniciales
# =============================================================================

groupID = '400138'
member = 'BCGsky'
    
sky_lim_pix = 720

   
# =============================================================================
# Parametros segun banda
# =============================================================================
bands = ['g', 'r', 'i']

color = {'g':'C0', 'r':'C2', 'i':'C3'}

SkyCorr = {'g':0,'r':0,'i':0}

Rlim = {'g':0,'r':0,'i':0}

subplot = {'g':131,'r':132,'i':133}

# =============================================================================
# Paths
# =============================================================================
# Path
masks_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/masks/'
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'
data_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/data/images/'+groupID


# =============================================================================
# Figura
# =============================================================================
plt.close('all')
plt.figure(1,figsize=(15,5))
plt.suptitle('Fitted sky')


# =============================================================================
# Loop sobre las bandas
# =============================================================================

for band in bands:
    
    path = results_path+'BCGsky/dist_prof_'+band+'-band.fits'
    
    prof = Table.read(path)

    # =============================================================================
    # Leemos todas las variables
    # =============================================================================
    
    r = prof['r']
    ct = prof['ct']
    ct_err = prof['ct err']

    maxsma = r[-1]


    # =============================================================================
    # Calculo del  cielo 
    # =============================================================================
    ind_lim = 700
    ind = np.array(r)>ind_lim
    #ind = ct-ct_err<=0+ct_err
    #ind = ind*np.isfinite(ct-ct_err)
    sky = np.nanmedian(ct[ind])
    SkyCorr[band] = sky


    
    # =============================================================================
    # Extraemos Rlim como el bin anterior a cuando se empieza a calcular el cielo
    # =============================================================================
    Rlim[band] = np.max(r[~ind])


    #==============================================================================
    # PLOT Perfil 
    #==============================================================================

    plt.subplot(subplot[band])

    # =============================================================================
    # Perfiles Ajustados
    # =============================================================================
    plt.title(band+'-band')
    # Perfil

    plt.errorbar(r, ct, yerr=ct_err, fmt='-o',markersize=4,zorder=0,color='C0')    
    plt.errorbar(r[ind], ct[ind], yerr=ct_err[ind], fmt='-o',markersize=4,zorder=1,color='C1')
    
    
    plt.axhline(0,color='k',lw=0.8)
    
    plt.axhline(sky,ls='dashed',color='C3',label = 'correction: '+str(sky))
    #plt.axvline(r[r>l[band]].min(),ls='dashed',color='r')
    
    plt.legend()
    

    plt.ylabel("Counts")
    plt.xlabel("$SMA$ [pix]")    
    plt.ylim(-0.0101,0.041)
    
    
    plt.tight_layout()
 
if member=='BCGsky':    
    plt.savefig(results_path+member+'/Dist_SkyCorr.pdf')
else:
    plt.savefig(results_path+member+'/Rlim.pdf')




# =============================================================================
# Buscamos el fichero de parametros y reemplazamos los de SkyCorr calculados
# =============================================================================

filename = results_path+'params.py'

if member=='BCGsky':    
    
    start = 'SkyCorr'
    replace = 'SkyCorr = {\'g\':'+str(SkyCorr['g'])+', \'r\':'+str(SkyCorr['r'])+', \'i\':'+str(SkyCorr['i'])+'}\n'
    
    import sys
    import fileinput
    for line in fileinput.input([filename], inplace=True):
        if line.strip().startswith(start):
            line = replace 
        sys.stdout.write(line)


start = 'Rlim_'+member
replace = 'Rlim_'+member+' = {\'g\':'+str(Rlim['g'])+', \'r\':'+str(Rlim['r'])+', \'i\':'+str(Rlim['i'])+'}\n'

import sys
import fileinput
for line in fileinput.input([filename], inplace=True):
    if line.strip().startswith(start):
        line = replace 
    sys.stdout.write(line)



# =============================================================================
# Restamos el cielo a la imagen 
#    ** para aplicar 2D models con imfit diretamente 
#    -- esto podria no ser definitivo
# ==============================================================================

for band in bands:
    data = fits.getdata(results_path+'PSF/PSFcorr-'+band+'.fits') 
    header = fits.getheader(results_path+'PSF/PSFcorr-'+band+'.fits') 
    dataSky = np.zeros_like(data)
    dataSky = data - SkyCorr[band]

    # =============================================================================
    # Guardamos la imagen con el cielo sustraido 
    # =============================================================================
    new_file = 'PSFcorr-'+band+'-SkySubDist.fits'
    save_fits_image(dataSky, header, results_path+'PSF/', new_file)


 
# =============================================================================
# PLOT Image of the distance-based apertures over the data
# =============================================================================
dist_map = fits.getdata(results_path+member+'/Dist_map.fits')

# Limit the values to show of the distances map to see the ellipses shape
distances_lim = np.copy(dist_map) # hay que cargar la imagen del mapa de distancias
distances_lim = distances_lim.astype("float")
distances_lim[distances_lim>sky_lim_pix] = np.nan

from matplotlib import colors
cmap = colors.ListedColormap(['k'])

plt.figure(2)

plt.imshow(data, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=5))
plt.imshow(distances_lim, cmap='tab20c', interpolation='nearest', origin='lower', alpha=0.7)
plt.savefig(results_path+member+'/Dist_Appertures.pdf')


   








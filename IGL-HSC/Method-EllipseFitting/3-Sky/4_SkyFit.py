import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"



# =============================================================================
# Parametros iniciales
# =============================================================================

groupID = '400138'
member = 'BCGsky'
    
    
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
    
    path = results_path+member+'/sbp_'+band+'-band.fits'
    
    sbp = Table.read(path)

    # =============================================================================
    # Leemos todas las variables
    # =============================================================================
    
    r = sbp['sma']
    sb = sbp['sb']
    err_up = sbp['err up']
    err_down = sbp['err down']

    
    ct = sbp['ct']
    ct_err = sbp['ct err']

    maxsma = r[-1]


    # =============================================================================
    # Calculo del  cielo 
    # =============================================================================
    ind = ct-ct_err<=0+ct_err
    ind = ind*np.isfinite(ct-ct_err)
    sky = np.nanmedian(ct[ind])
    SkyCorr[band] = sky

    
    
    # =============================================================================
    # Extraemos Rlim como el bin anterior a cuando se empieza a calcular el cielo
    # =============================================================================
    Rlim[band] = np.max(r[~ind])


    #==============================================================================
    # PLOT Perfil 
    #==============================================================================


    # =============================================================================
    # Pintamos los datos
    # =============================================================================

    plt.subplot(subplot[band])

    # =============================================================================
    # Perfiles Ajustados
    # =============================================================================
    plt.title(band+'-band')
    # Perfil

    plt.errorbar(r, ct, yerr=ct_err, fmt='-o',markersize=4,zorder=0,color='C0')    
    plt.errorbar(r[ind], ct[ind], yerr=ct_err[ind], fmt='-o',markersize=4,zorder=1,color='C1')
    
    
    plt.axhline(0,color='k',lw=0.8)
    
    if member=='BCGsky':  
        plt.axhline(sky,ls='dashed',color='C3',label = 'correction: '+str(sky))
    #plt.axvline(r[r>l[band]].min(),ls='dashed',color='r')
    
    plt.legend()
    

    plt.ylabel("Counts")
    plt.xlabel("$sma$ [pix]")    
    plt.ylim(-0.0101,0.041)
    
    
    plt.tight_layout()
 
if member=='BCGsky':    
    plt.savefig(results_path+member+'/SkyCorr.pdf')
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
# =============================================================================

#SkyCorr = {'g':-0.000364439, 'r':-0.0016048, 'i':-0.000529909}

for band in bands:
    #data = fits.getdata(data_path+'-'+band+'.fits')
    #header = fits.getheader(data_path+'-'+band+'.fits')
    data = fits.getdata(results_path+'PSF/PSFcorr-'+band+'.fits') 
    header = fits.getheader(results_path+'PSF/PSFcorr-'+band+'.fits') 
    dataSky = np.zeros_like(data)
    dataSky = data - SkyCorr[band]

    # =============================================================================
    # Guardamos la imagen con el cielo sustraido 
    # =============================================================================
    
    masked_hdu = fits.PrimaryHDU(dataSky)  # we create a PrimaryHDU object to encapsulate the data:
    wcs = WCS(header)
    masked_hdu.header.update(wcs.to_header())
    masked_hdu.header.append(('HISTORY', 'Sky subtracted', 'C. Martinez-Lombilla'), end=True)
    
    # !!!!! -- cambiar nombre de la imagen??
    masked_hdu.writeto(results_path+'PSF/PSFcorr-'+band+'-SkySub.fits', overwrite=True) # write to a new file
    #masked_hdu.writeto(results_path+member+'/'+band+'-SkySub.fits', overwrite=True) # write to a new file




import os, glob, sys, warnings, array, re, math, time, copy
import numpy                  as np
import matplotlib.pyplot      as plt
import matplotlib.cm          as cm
from   astropy.io             import fits, ascii
from   astropy.wcs            import WCS
from   astropy                import units as u
from   astropy.coordinates    import SkyCoord
from astropy.stats            import sigma_clipped_stats
#######

mag_isophot = np.arange(20, 27, 0.2)
groupID = '400138'
groupID = '400138'
members_core = ['BCG', '1660730', '1660646', '1660615']
bands = ['g', 'r', 'i']

# =============================================================================
# Paths
# =============================================================================
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'


# =============================================================================
# Fichero de parametros
# =============================================================================
# Copiamos el fichero de parametros para asegurarnos de que es el mismo

cmd = 'cp '+results_path+'params.py .'

os.system(cmd)
import params

# Parametros del grupo
z = params.z

# Parametros del instrumento y la imagen
zp = params.zp
pixscale = params.pixelscale
SkyRMS = params.SkyRMS


# A y Correccion K
k_corrs = params.k_corrs
extinctions = params.extinctions
SkyRMS = params.SkyRMS

# sb dimming: -2.5*np.log10(intens) - 10*np.log10(1+z) -- pq el dimming va como (1+z)**-4 -- Calvi+2014  doi: 10.1088/0004-637X/796/2/102
dimming = 10*np.log10(1+z) 

#Posiciones centrales del core
Xcore = []
Ycore = []
for memberc in members_core:
    Xcore.append(int(round(params.Xc[memberc],0)))
    Ycore.append(int(round(params.Yc[memberc],0)))

    
os.system('rm params.py')





# Images
image_name = ['CoreMasked-g.fits', 'CoreMasked-r.fits', 'CoreMasked-i.fits' ]


## For the reference image :

im, hdr = fits.getdata(results_path+'masks/'+image_name[0], header = True)

isophot_counts = []
test           = np.zeros_like(im)

for tt in np.arange(0, len(image_name)):
    counts = 10**(-(mag_isophot- zp - 5.*np.log10(pixscale) + dimming + k_corrs[bands[tt]] + extinctions[bands[tt]])/2.5)
    image, hdr = fits.getdata(results_path+'masks/'+image_name[tt], header = True)
    result      = []

    for pp in np.arange(0, len(mag_isophot)-1):
    ## Define la isofota en banda X
        ind_isophot = (im < counts[pp]) & (im > counts[pp+1])

        #Test
        test[ind_isophot] = pp + tt
        mean_counts = sigma_clipped_stats(image[ind_isophot], sigma = 3)[0]
        result.append(mean_counts)

    isophot_counts.append(result)
    fits.writeto(results_path+'tests/'+'isophotes_'+ str(tt) + '.fits', test, overwrite = True)


sys.exit()

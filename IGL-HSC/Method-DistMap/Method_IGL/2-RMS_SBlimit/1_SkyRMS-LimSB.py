# Basics
import os
import numpy as np


import matplotlib.pyplot as plt
# Astropy
from astropy.io import fits

import sys
sys.path.append('/Users/felipe/Dropbox/CML/') # path CML.py Feli
sys.path.append('/Users/z3530031/Dropbox/CML/') # path CML.py Cris

sys.path.append('/Users/felipe/Dropbox/CML/WHT/') # path del fichero de parametros Feli
sys.path.append('/Users/z3530031/Dropbox/CML/WHT/') # path del fichero de parametros Cris

# CML
from CML import RectStat

#==============================================================================
# Parametros iniciales
#==============================================================================
groupID = '400138'


# =============================================================================
# Diccionarios de las bandas
# =============================================================================

bands = ['g','r','i'] 

# =============================================================================
# Paths
# =============================================================================
masks_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/masks/'
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'



# =============================================================================
# Cargamos los parametros del fichero de parametros
# =============================================================================

cmd = 'cp '+results_path+'params.py .'

os.system(cmd)
import params


# Parametros del instrumento
zp = params.zp
pixelscale = params.pixelscale

os.system('rm params.py')


# =============================================================================
# Tamano de la apertura
# =============================================================================
# En arcsec
size = 10 

# En pixels
size_pix = size/pixelscale

# Numero de aperturas aleatorias
n = 2500


# =============================================================================
# Bucle sobre las bandas
# =============================================================================
SkyRMS = {'g': 0, 'r': 0, 'i': 0}
sb_limit = {'g':0, 'r':0, 'i':0}

for band in bands:
    # =============================================================================
    # Cargamos la imagen, la mascara y la enmascaramos
    # =============================================================================
    data = fits.getdata(results_path+'PSF/PSFcorr-'+band+'.fits')
    mask = fits.getdata(masks_path+'RMSMask-'+band+'.fits')
    data[np.isnan(mask)] = np.nan 

    # Vamos a definir un cuadrado del tamano de la imagen centrada en el centro de la misma
    w = size_pix
    h = size_pix

    # =============================================================================
    # Calculo del rms
    # =============================================================================
    
    r = RectStat(data, w, h, n=n)
        
    mean, median, stddev = r.stat()    # rms = stddev obtenidos con sigma_clipped_stats
    
    plt.figure(bands.index(band))
    r.hist(bins=100)
    
   
    # =============================================================================
    # Definicion Tesis Rulo  +  Roman+20
    # =============================================================================    

    mu = -2.5 * np.log10(3*stddev/(pixelscale * size)) + zp #Roman+20(galactic cirri, Appendice A)

    SkyRMS[band] = stddev
    sb_limit[band] = np.round(mu,2)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           



# =============================================================================
# Buscamos el fichero de parametros y reemplazamos los de sb_limit y SkyRMS calculados
# =============================================================================
filename = results_path+'params.py'

#SkyRMS
start = 'SkyRMS'
replace = 'SkyRMS = {\'g\':'+str(SkyRMS['g'])+', \'r\':'+str(SkyRMS['r'])+', \'i\':'+str(SkyRMS['i'])+'}\n'

import sys
import fileinput
for line in fileinput.input([filename], inplace=True):
    if line.strip().startswith(start):
        line = replace 
    sys.stdout.write(line)


# sl_limit
start = 'sb_limit'
replace = 'sb_limit = {\'g\':'+str(sb_limit['g'])+', \'r\':'+str(sb_limit['r'])+', \'i\':'+str(sb_limit['i'])+'}\n'

for line in fileinput.input([filename], inplace=True):
    if line.strip().startswith(start):
        line = replace 
    sys.stdout.write(line)


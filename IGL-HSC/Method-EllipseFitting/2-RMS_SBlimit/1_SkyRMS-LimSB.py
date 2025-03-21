# Basics
import os
import numpy as np


import matplotlib.pyplot as plt
# Astropy
from astropy.io import fits

import fileinput
import sys

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
sb_limit_3s = {'g':0, 'r':0, 'i':0}
sb_limit_1s = {'g':0, 'r':0, 'i':0}


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
        
    mean, median, stddev = r.stat()    # rms = stddev obtained using sigma_clipped_stats
    
    plt.figure(bands.index(band))  # to check that the sky pixel values are gaussian
    r.hist(bins=100)
    
   
    # =============================================================================
    # Definition from Raul's thesis and Javi (Roman+20)
    # =============================================================================    

    mu_s3 = -2.5 * np.log10(3*stddev/(pixelscale * size)) + zp #Roman+20(galactic cirri, Appendice A)
    mu_s1 = -2.5 * np.log10(1*stddev/(pixelscale * size)) + zp #Roman+20(galactic cirri, Appendice A)

    SkyRMS[band] = stddev
    
    sb_limit_3s[band] = np.round(mu_s3,2)
    sb_limit_1s[band] = np.round(mu_s1,2)



# =============================================================================
# Buscamos el fichero de parametros y reemplazamos los de sb_limit y SkyRMS calculados
# =============================================================================
filename = results_path+'params.py'

#SkyRMS
start = 'SkyRMS'
replace = 'SkyRMS = {\'g\':'+str(SkyRMS['g'])+', \'r\':'+str(SkyRMS['r'])+', \'i\':'+str(SkyRMS['i'])+'}\n'


for line in fileinput.input([filename], inplace=True):
    if line.strip().startswith(start):
        line = replace 
    sys.stdout.write(line)


# sb_limit
start = 'sb_limit'
replace = 'sb_limit = {\'g\':'+str(sb_limit_3s['g'])+', \'r\':'+str(sb_limit_3s['r'])+', \'i\':'+str(sb_limit_3s['i'])+'}\n'

for line in fileinput.input([filename], inplace=True):
    if line.strip().startswith(start):
        line = replace 
    sys.stdout.write(line)


# sb_limit 1 sigma
start = 'sb_limit_1sigma'
replace = 'sb_limit_1sigma = {\'g\':'+str(sb_limit_1s['g'])+', \'r\':'+str(sb_limit_1s['r'])+', \'i\':'+str(sb_limit_1s['i'])+'}\n'

for line in fileinput.input([filename], inplace=True):
    if line.strip().startswith(start):
        line = replace 
    sys.stdout.write(line)


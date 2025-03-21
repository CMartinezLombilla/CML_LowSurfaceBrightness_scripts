import pickle
import numpy as np
import numpy.ma as ma
import os
from astropy.io import fits
from astropy.table import Table

# =============================================================================
# importamos CML
# =============================================================================

import sys
sys.path.append('/Users/felipe/Dropbox/CML/') # path del CML Feli
sys.path.append('/Users/z3530031/Dropbox/CML/') # path del CML Cris

sys.path.append('/Users/felipe/Dropbox/CML/WHT/') # path del fichero de parametros Feli
sys.path.append('/Users/z3530031/Dropbox/CML/WHT/') # path del fichero de parametros Cris



#==============================================================================
#Initial parameters
#==============================================================================

GroupID = '400138'
member = 'BCGsky'


# =============================================================================
# Diccionario de las banda
# =============================================================================
bands = ['g', 'r', 'i']


# =============================================================================
# Paths
# =============================================================================
# Path
masks_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+GroupID+'/masks/'
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+GroupID+'/'
data_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/data/images/'+GroupID


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
pixelscale = params.pixelscale

# Parametros por banda
SkyCorr = params.SkyCorr
SkyRMS = params.SkyRMS
extinctions = params.extinctions
k_corrs = params.k_corrs
# sb dimming: -2.5*np.log10(intens) - 10*np.log10(1+z) -- pq el dimming va como (1+z)**-4 -- Calvi+2014  doi: 10.1088/0004-637X/796/2/102
dimming = 10*np.log10(1+z) 

# Parametros por galaxia
mask_file = params.mask_file[member]
os.system('rm params.py')



# =============================================================================
# Print check
# =============================================================================
print('# =================================================')
print('               Member: ' +str(member)+'               ')
print('# =================================================')



# =============================================================================
# Cargamos el objeto que continene el ajuste original
# =============================================================================
elipses = pickle.load(open(results_path+member+'/'+"main_model.p", "rb" ) )


# =============================================================================
# Loop sobre todas las bandas
# =============================================================================
for band in bands:
    # =============================================================================
    # Leemos los datos y los enmascaramo 
    # =============================================================================
    
    if member == 'BCGsky': # No aplicamos la correcion si lo estamos haciendo para el cielo 
        SkyCorr[band] = 0 
    
    mask = fits.getdata(masks_path+mask_file+'-'+band+'.fits')       
    #data = fits.getdata(data_path+'-'+band+'.fits')
    data = fits.getdata(results_path+'PSF/PSFcorr-'+band+'.fits')-SkyCorr[band]
    data_masked  = ma.masked_array(data, mask)
    
    k_corr = k_corrs[band]
    A = extinctions[band]
    
    sbp = elipses.run_over(data_masked,zp,pixelscale, rms_sky=SkyRMS[band], A=A, dimming=dimming, k_corr=k_corr)


    sbp.write(results_path+member+'/sbp_'+band+'-band.fits', overwrite=True)
        
"""
index = sbp['SNR']>0.001
plt.close('all')
plt.scatter(sbp['sma'][index],sbp['sb'][index],label='all')
plt.scatter(sbp['sma'][index],- 2.5*np.log10(sbp['ct'][index]) + zp ,label='zp')
plt.scatter(sbp['sma'][index],- 2.5*np.log10(sbp['ct'][index]) + zp - A,label='zp, A')
plt.scatter(sbp['sma'][index],- 2.5*np.log10(sbp['ct'][index]) + zp - dimming,label='zp, dimming' )
plt.scatter(sbp['sma'][index],- 2.5*np.log10(sbp['ct'][index]) + zp - k_corr,label='zp, k_coor')
#plt.axhline(27.5,ls='dashed',color='k')
plt.gca().invert_yaxis()
plt.ylabel('$\mu_{'+band+'}$ (mag/arsec$^{2}$)',size=14)
plt.title(band+'-band', y=0.9,x=0.5,size=15)
plt.xlabel('SMA (pix)',size=14)
plt.legend()
plt.show()

"""




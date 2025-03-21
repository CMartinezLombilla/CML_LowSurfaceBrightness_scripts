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
from regions import read_ds9
from astropy.stats import sigma_clipped_stats


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
from CML import save_fits_image
from CML import Counts2Sb



# =============================================================================
# Parametros iniciales del grupo/imagen/IGL/aperturas
# =============================================================================
groupID = '400138'
bands = ['g', 'r', 'i']


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


os.system('rm params.py')

# Images
im_names = ['CoreMasked-g.fits', 'CoreMasked-r.fits', 'CoreMasked-i.fits' ]
im, hdr = fits.getdata(results_path+'masks/'+im_names[1], header=True)
wcs = WCS(hdr)
shape = im.shape


# Region for the IGL core:
core_reg = read_ds9(results_path+'IGL/StellPop/ds9Regs/IGL_core.reg')
mask_core_reg = np.zeros_like(im)

for i in core_reg:    
   mask_core_reg = mask_core_reg +i.to_pixel(wcs).to_mask('center').to_image(shape)
mask_core_reg[mask_core_reg==0] = np.nan


ct_median_core = []
for im_name in im_names: 
    im = fits.getdata(results_path+'masks/'+im_name)
    ct_median_core.append(sigma_clipped_stats(im[~np.isnan(mask_core_reg)], sigma = 1)[1])



# Region for the IGL tail:
tail_reg = read_ds9(results_path+'IGL/StellPop/ds9Regs/IGL_tail.reg')
mask_tail_reg = np.zeros_like(im)

for i in tail_reg:    
   mask_tail_reg = mask_tail_reg +i.to_pixel(wcs).to_mask('center').to_image(shape)
mask_tail_reg[mask_tail_reg==0] = np.nan


ct_median_tail = []
for im_name in im_names: 
    im = fits.getdata(results_path+'masks/'+im_name)
    ct_median_tail.append(sigma_clipped_stats(im[~np.isnan(mask_tail_reg)], sigma = 1)[1])



# INPUT DATA: structures and counts to convert
structures = ['BCG', '1660545', '1660615', '1660646', '1660730', '2301069', 'IGL_core', 'IGL_tail']

ct_g = np.array([4.14E-01, 1.63E-01, 8.83E-02, 3.86E-01, 2.78E-01, 0.0880229, ct_median_core[0], ct_median_tail[0]])
ct_r = np.array([1.2151, 0.406928, 0.246497, 1.16597, 0.805679, 0.209255, ct_median_core[1], ct_median_tail[1]])
ct_i = np.array([1.98156, 0.613439, 0.400243, 1.89597, 1.27768, 0.368664, ct_median_core[2], ct_median_tail[2]])



# Obtain the SB values corrected by k-corr, gal ext. and dimming
c2sb_g = Counts2Sb(zp, pixscale, A=extinctions['g'], dimming=dimming, k_corr=k_corrs['g'])
sb_g = c2sb_g.do(ct_g)

c2sb_r = Counts2Sb(zp, pixscale, A=extinctions['r'], dimming=dimming, k_corr=k_corrs['r'])
sb_r = c2sb_r.do(ct_r)

c2sb_i = Counts2Sb(zp, pixscale, A=extinctions['i'], dimming=dimming, k_corr=k_corrs['i'])
sb_i = c2sb_i.do(ct_i)


# RESULTS: COLORS corrected
my_g_r = sb_g - sb_r
my_r_i = sb_r - sb_i












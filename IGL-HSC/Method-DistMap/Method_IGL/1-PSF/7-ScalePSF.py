import numpy as np
import numpy.ma as ma
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.close('all')
plt.ioff()

plt.rcParams["font.family"] = "serif"
plt.rc('text', usetex=False)

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

from regions import read_ds9

# =============================================================================
# importamos CML
# =============================================================================

import sys
sys.path.append('/Users/felipe/Dropbox/CML/') # path del CML Feli
sys.path.append('/Users/z3530031/Dropbox/CML/') # path del CML Cris

sys.path.append('/Users/felipe/Dropbox/CML/WHT/') # path del fichero de parametros Feli
sys.path.append('/Users/z3530031/Dropbox/CML/WHT/') # path del fichero de parametros Cris

from CML import LogNorm
from CML import ScalePSF


# =============================================================================
# Fichero de parametros
# =============================================================================
# Copiamos el fichero de parametros para asegurarnos de que es el mismo

groupID = '400138'
data_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/data/'
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'
    
cmd = 'cp '+results_path+'params.py .'

os.system(cmd)
import params

# Parametros del grupo
bcgID = params.bcgID
z = params.z

SkyPSF = params.SkyPSF
    
os.system('rm params.py')

# =============================================================================
# Parameters
# =============================================================================

bands = ['g']#,'r','i']

f_min = {'g':abs(SkyPSF['g']*2), 
         'r':abs(SkyPSF['r']*2), 
         'i':abs(SkyPSF['i']*2)} # n sigmas of the SkyRMS

f_max= {'g':0.04, 'r':0.04, 'i':0.04} # Sat limit * 0.01

n_bsiter = 2000

scale_factor={'g':{1:285077, 2:50000, 3:0}, 
              'r':{1:335333, 2:70790, 3:0}, 
              'i':{1:604137, 2:76562, 3:0}}
#'g':{1:275077, 2:50000, 3:0}

for band in bands:
    stars = [1,2,3]
    
    star_center = {1: [1599.4475,650.55127],
                   2: [1732.4475,772.55127],
                   3:[306.52006,1738.973]}
    
    
    model = 0
    
    # =============================================================================
    # =============================================================================
    
    
    # Data
    load_data = fits.open(data_path+'images/'+groupID+'-'+band+'.fits')
    #load_data = fits.open('/Users/felipe/OneDrive - UNSW/Analysis-HSC/data/images/400138-'+band+'.fits')
    data = load_data[1].data
    header = load_data[1].header
    wcs_data = WCS(header)

    
    
    
    # Hacemos un for sobre las estrellas
    
    for star in stars:
        # Mask
        mask = fits.getdata(results_path+'PSF/masks/mask_s'+str(star)+'_final.fits')
        data_masked  = ma.masked_array(data - model, mask)
        
        # PSF
        psf = fits.getdata(data_path+'brightStars/PSFsBaena/psf_outer-inter-inner_'+band.capitalize()+'.fits')
        psf_hed = fits.getheader(data_path+'brightStars/PSFsBaena/psf_outer-inter-inner_'+band.capitalize()+'.fits')
        wcs = WCS(psf_hed)
        regions = read_ds9(data_path+'brightStars/PSFsBaena/ds9Regions/PSFspikes-'+band+'.reg')
        shape = psf.shape
        mask = np.zeros_like(psf)
        
        for ii in regions:    
            mask = mask +ii.to_mask('center').to_image(shape)
            
        
        # =============================================================================
        # Escalado de la PSF y modelo 2D
        # =============================================================================
        scale = ScalePSF(data_masked, psf, star_center=star_center[star], psf_mask=mask, f_min=f_min[band], f_max = f_max[band], r_min=25, r_max=450, step=20)
        
        #model = model + scale.model_optimize(scale_factor=0,n_bs=n_bsiter)
        
        model = model + scale.model_optimize(scale_factor=scale_factor[band][star],n_bs=n_bsiter)
        
        # =============================================================================
        #  Graficas de control
        # =============================================================================
        
        
        plt.close('all')
         
        scale.plot(savefig=results_path+'PSF/'+band+'-band_star'+str(star)+'.pdf')
            
    
    # =============================================================================
    # Guardamos el del scaterred light field con todas las estrellas a sustraer
    # =============================================================================
    model_hdu = fits.PrimaryHDU(model)  # we create a PrimaryHDU object to encapsulate the data:
    model_hdu.header.update(wcs_data.to_header())
    model_hdu.writeto(results_path+'PSF/scatter_field-'+band+'.fits', overwrite=True) # write to a new file    
        
        
    plt.figure(2)
    plt.imshow(data - model, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=10))
    plt.savefig(results_path+'PSF/'+band+'-band_res.pdf')    



    # =============================================================================
    # Guardamos la imagen con las estrellas brillantes sustraidas
    # =============================================================================
    im_subpsf = data - model
    model_hdu = fits.PrimaryHDU(im_subpsf)  # we create a PrimaryHDU object to encapsulate the data:
    model_hdu.header.update(wcs_data.to_header())
    model_hdu.writeto(results_path+'PSF/PSFcorr-'+band+'.fits', overwrite=True) # write to a new file    
        





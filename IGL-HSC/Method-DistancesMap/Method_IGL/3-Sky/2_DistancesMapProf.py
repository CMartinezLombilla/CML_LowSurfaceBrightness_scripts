import numpy as np
import os, sys
from astropy.io import fits
from astropy.wcs import WCS

import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.close('all')

# =============================================================================
# import CML
# =============================================================================

import sys
sys.path.append('/Users/z3530031/Dropbox/CML/') # path del CML Cris

from CML import LogNorm
from CML import Pix2Kpc
from CML import Kpc2Pix
from CML import dist_ellipse_map
from CML import dist_ellipse_ct_prof
from CML import save_fits_image


#==============================================================================
# Initial params -- EDIT AS REQUIRED
#==============================================================================
GroupID = '400138'
members = ['BCG', '1660646','1660730',] # Galaxies I want to get distances from

# Bands
bands = ['g', 'r', 'i']   

# Params profile
sma0 = 2
nbins = 27


# =============================================================================
# Main path for images and masks -- EDIT AS REQUIRED
# =============================================================================
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+GroupID+'/'
data_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/data/images/'+GroupID


# =============================================================================
# Params file
# =============================================================================
# We copy the file to make sure is always the same

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
dimming = 10*np.log10(1+z) 


# Sky rms and mask
SkyRMS = params.SkyRMS
mask_file = params.mask_file['BCG']

# Params and lists of params for each member
maxsma = params.maxsma['BCGsky']
step = params.step['BCGsky']

x_gals = []
y_gals = []
ell_gals = []
pa_gals = [] # PA 2d-imfit + 90 deg


for member in members:              
    # Galaxy params
    x_gals.append(params.Xc[member])
    y_gals.append(params.Yc[member])
    ell_gals.append(params.ell[member])
    pa_gals.append(params.theta[member]) # deg. (Desde el eje x, counterclockwise) esto sería nuestro theta


os.system('rm params.py')



# =============================================================================
# Code to get the distances map
# ============================================================================= 

# Load an image in any band (just want the size and header... same at all bands!)
im, hdr = fits.getdata(results_path+'PSF/PSFcorr-g.fits', header = True)

# Get the distances map
ny,nx = np.shape(im)
dist = []

for tt in np.arange(len(members)):
        grid = dist_ellipse_map([nx,ny], x_gals[tt], y_gals[tt], 1-ell_gals[tt], pa_gals[tt])
        grid = np.array(grid, dtype = int)
        dist.append(grid)

dist_map = np.amin(dist, axis = 0) #* 0.06 * arcsec_to_kpc[ii]
  
# Check distance map plot
plt.imshow(np.log10(dist_map), origin='lower'), plt.show()

# Save distances map .fits
new_file = 'Dist_map.fits'
save_fits_image(dist_map, hdr, results_path+'BCGsky/', new_file)


  
# =============================================================================
# Extract and save the (counts) profile using the distances map
# =============================================================================

for band in bands:
          
    # Load the images and masks at each band
    mask = fits.getdata(results_path+'masks/'+mask_file+'-'+band+'.fits')   
    data, hdr = fits.getdata(results_path+'PSF/PSFcorr-'+band+'.fits', header = True)
    
    # Masking the data    
    im  = np.copy(data)
    im[np.isnan(mask)] = np.nan
    

    ct_prof = dist_ellipse_ct_prof(im, dist_map, sma0, step, nbins, SkyRMS[band], zp, pixscale, A=extinctions[band], dimming=dimming, k_corr=k_corrs[band])
    ct_prof.write(results_path+'BCGsky/dist_prof_'+band+'-band.fits', overwrite=True)




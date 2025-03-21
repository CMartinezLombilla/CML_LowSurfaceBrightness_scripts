import numpy as np
import os, sys, glob, copy
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats        

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



#==============================================================================
# Initial params -- EDIT AS REQUIRED
#==============================================================================
GroupID = '400138'
members = ['BCG', '1660646','1660730',] # Galaxies I want to get distances from

# Fix ellipse centrer?
fix_center = True


# Bands
bands = ['g', 'r', 'i']   

# Plots 
subplot = {'g':131,'r':132,'i':133}


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
pixelscale = params.pixelscale

# SB Profile Corrections
k_corrs = params.k_corrs
extinctions = params.extinctions
# sb dimming: -2.5*np.log10(intens) - 10*np.log10(1+z) -- Calvi+2014  doi: 10.1088/0004-637X/796/2/102
dimming = 10*np.log10(1+z) 
mask_file = params.mask_file['BCG']

os.system('rm params.py')



# =============================================================================
# Function: returns distances map from the centre of each galaxy [Mireia]
# =============================================================================


def dist_ellipse(N, xc, yc, ratio, pos_ang):
    import numpy as np
    #Convert POS_ANG to radians
    ang = pos_ang * np.pi / 180.
    cosang = np.cos(ang)
    sinang = np.sin(ang)
    if len(N) == 2 :
        nx = N[0]
        ny = N[1]
    elif len(N)==1 :
        ny = N[0]
        nx = N[0]
    x = np.arange(0, nx,1) - [xc-1.]*nx
    y = np.arange(0, ny,1) - [yc-1.]*ny
    im = np.zeros([ny, nx])
    #                     ;Rotate pixels to match ellipse orientation
    xcosang = x*cosang
    xsinang = x*sinang
    for ii in np.arange(0, ny, 1):
         xtemp =  xcosang + y[ii]*sinang
         ytemp = -xsinang + y[ii]*cosang
         im[ii,:] = np.sqrt((ytemp/ratio)**2 + xtemp**2 )
    return im


# =============================================================================
# Code to get the distances map
# =============================================================================

#member = 'BCG'           
cmd = 'cp '+results_path+'params.py .'
os.system(cmd)
import params

x_gals = []
y_gals = []
ratio_gals = []
pa_gals = [] # PA 2d-imfit + 90 deg


for member in members:           
    
    # Galaxy params
    x_gals.append(params.Xc[member])
    y_gals.append(params.Yc[member])
    ratio_gals.append(params.ell[member])
    pa_gals.append(params.theta[member]) # deg. (Desde el eje x, counterclockwise) esto sería nuestro theta

os.system('rm params.py')



# Some plot set ups...
plt.close('all')
plt.figure(1,figsize=(15,5))
plt.suptitle('Fitted sky')



# Loop over the bands
for band in bands:
          
    # Load the images and masks
    mask = fits.getdata(results_path+'masks/'+mask_file+'-'+band+'.fits')   
    data, hdr = fits.getdata(results_path+'PSF/PSFcorr-'+band+'.fits', header = True)
    
    # Masking the data    
    im  = np.copy(data)
    im[np.isnan(mask)] = np.nan
    
    
    ny,nx = np.shape(im)
    dist = []
    
    for tt in np.arange(len(members)):
            grid = dist_ellipse([nx,ny], x_gals[tt], y_gals[tt], 1-ratio_gals[tt], pa_gals[tt])
            grid = np.array(grid, dtype = int)
            dist.append(grid)
    
    distances = np.amin(dist, axis = 0) #* 0.06 * arcsec_to_kpc[ii]
    
    
    
    
    # Check distance map plot
    #plt.imshow(np.log10(distances), origin='lower'), plt.show()
    
    
    
    # =============================================================================
    # Code to extract the (counts) profile using the distances map
    # =============================================================================
     
    maxsma = params.maxsma['BCGsky']
    step = 0.25
    sma0 = 2
    nbins = 25
    
    ell_bin = [sma0*(1. + step)]
    ct = []
    r = []
    
    for ii in range(nbins):
        ell_bin.append(ell_bin[ii]*(1. + step))
        im_ell = im[(distances >= ell_bin[ii]) & (distances < ell_bin[ii+1])]
        ct.append(sigma_clipped_stats(im_ell, sigma=3, mask_value=np.nan)[1])
        r.append(ell_bin[ii] + (ell_bin[ii]+ ell_bin[ii+1])/2)
    
    
    # =============================================================================
    # Sky determination
    # =============================================================================
    ind_lim = 700
    ind = np.array(r)>ind_lim

    sky = np.nanmedian(np.array(ct)[ind])
    SkyCorr = sky
    
    
    
    #==============================================================================
    # PLOT profiles 
    #==============================================================================
    plt.subplot(subplot[band])
    
    # =============================================================================
    # Perfiles Ajustados
    # =============================================================================
    plt.title(band+'-band')
    
    
    plt.plot(r, ct, marker='o', markersize=4,zorder=1,color='C0') 
    plt.plot(np.array(r)[ind], np.array(ct)[ind], linestyle='dashed', marker='o', markersize=4,zorder=1,color='C1')
    
    
    plt.axhline(0,color='k',lw=0.8)
    
    plt.axhline(sky,ls='dashed',color='C3',label = 'correction: '+ str(round(sky,6)))
    
    plt.legend()
    
    
    plt.ylabel("Counts")
    plt.xlabel("Semi-minor axis [pix]")    
    plt.ylim(-0.0101,0.041)
    
    
plt.tight_layout()
    
    

# Limit the values to show of the distances image to see the ellipses shape
distances_lim = np.copy(distances)
distances_lim = distances_lim.astype("float")
distances_lim[distances_lim>ind_lim] = np.nan

from matplotlib import colors
cmap = colors.ListedColormap(['k'])

plt.figure(2)

plt.imshow(im, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=5))
plt.imshow(distances_lim, cmap='tab20c', interpolation='nearest', origin='lower', alpha=0.7)












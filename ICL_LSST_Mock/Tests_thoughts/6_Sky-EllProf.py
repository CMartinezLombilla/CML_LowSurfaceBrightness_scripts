import numpy as np
import numpy.ma as ma
import os
import glob
import matplotlib.pyplot as plt
plt.close('all')
plt.ioff()

plt.rcParams["font.family"] = "serif"
plt.rc('text', usetex=False)

from astropy.io import fits



# =============================================================================
# importamos CML
# =============================================================================

from CML import EllipticalProfile
from CML import LogNorm
from CML import dist_ellipse_map
from CML import save_fits_image
from CML import dist_ellipse_ct_prof

# Division in "base" whatever
def baseround(x, base=10):
    return base * round(x/base)




# =============================================================================
# Initial params
# =============================================================================
sim = 'Horizon_AGN'
cluster = '048'
orientation = 'xy'

pixscale = 0.4
zp = 0

# Fix ellipse centrer?
fix_center = True


# =============================================================================
# Paths 
# =============================================================================
results_path = '/Users/z3530031/unsw/Work/Projects/LSST/Chall3-ICL/results/'+sim+'/'+cluster+'_'+orientation+'/'
mask_path = results_path+'Masks/'


# =============================================================================
# Create a folder to save script products
# =============================================================================
try:
    # Create target Directory
    os.mkdir(results_path+'Sky')
except:
    pass



# =============================================================================
# Load images and masks -- EDIT AS REQUIRED
# =============================================================================

data, hdr = fits.getdata(results_path+'Images/00761_0000048_0.05_xy_r_4Mpc_bin2_x1e+12.fits', header = True) 
#data, hdr = fits.getdata(mask_path+'final-'+sim+'-'+cluster+'-'+orientation+'.fits', header = True) 

mask = np.zeros_like(data)
mask[np.isnan(data)] = np.nan


# =============================================================================
# BCG Params for ellipse fitting
# =============================================================================

size = np.size(data,1)
Xc = size/2
Yc = size/2
ell = 0.3
pa = 1
sma0 = 15
step = 0.2 # default is 0.1
nbins = 30


# =============================================================================
# Code to get the distances map
# ============================================================================= 

# Get the distances map
ny,nx = np.shape(data)
dist = []


grid = dist_ellipse_map([nx,ny], Xc, Yc, 1-ell, pa)
grid = np.array(grid, dtype = int)
dist.append(grid)

dist_map = np.amin(dist, axis = 0) #* 0.06 * arcsec_to_kpc[ii]
  
# Check distance map plot
plt.imshow(np.log10(dist_map), origin='lower'), plt.show()

# Save distances map .fits
new_file = 'Dist_map.fits'
save_fits_image(dist_map, hdr, results_path+'Sky/', new_file)


  
# =============================================================================
# Extract and save the (counts) profile using the distances map
# =============================================================================

ct_prof = dist_ellipse_ct_prof(data, dist_map, sma0, step, nbins, 0, zp, pixscale, A=0, dimming=0, k_corr=0)
ct_prof.write(results_path+'Sky/sky_dist_prof.fits', overwrite=True)


# =============================================================================
# Leemos todas las variables
# =============================================================================

r = ct_prof['r']
ct = ct_prof['ct']
ct_err = ct_prof['ct err']

maxsma = r[-1]


# =============================================================================
# Calculo del  cielo 
# =============================================================================
ind_lim = 6000
ind = np.array(r)>ind_lim

sky = np.nanmedian(ct[ind])




#==============================================================================
# PLOT Perfil 
#==============================================================================

plt.figure(2)

# =============================================================================
# Perfiles Ajustados
# =============================================================================
plt.title(sim+'-'+cluster+'-'+orientation)
# Perfil

plt.errorbar(r, ct, yerr=ct_err, fmt='-o',markersize=4,zorder=0,color='C0')    
plt.axhline(0,color='k',lw=0.8)
plt.axhline(sky,ls='dashed',color='C3',label = 'correction: '+str(sky))
#plt.axvline(r[r>l[band]].min(),ls='dashed',color='r')

plt.legend()


plt.ylabel("Counts")
plt.xlabel("$SMA$ [pix]")    
plt.ylim(-0.01,0.1)

plt.tight_layout()
 
plt.savefig(results_path+'Sky/SkyCorr.pdf')





'''
# =============================================================================
# Restamos el cielo a las imagenes
# ==============================================================================
sims = ['Horizon_AGN', 'Hydrangea', 'Magneticum', 'TNG_100'] 
orientations = ['xy', 'yz', 'xz']
bands = ['r']

npix_bin = 2 #number of pixels for the binning
factor = 1e12 # scale factor for multiplying the data values
factor_str = format(factor, '.0e')


for sim in sims:
    
    # =============================================================================
    # Paths 
    # =============================================================================
    data_dir = '/Users/z3530031/unsw/Work/Projects/LSST/Chall3-ICL/data/'+sim+'/'
    res_dir = '/Users/z3530031/unsw/Work/Projects/LSST/Chall3-ICL/results/'
    
    #image_file = '00761_0000048_0.05_xy_r_4Mpc'
    
    image_files = glob.glob(res_dir+'*_r_4Mpc.fits') # Load a list 
    
    for clust_files in image_files:
        
        cluster = clust_files[-23:-20]
        orientation = clust_files[-14:-12]
        

        #res_dir+sim+'/'+cluster+'_'+orientation+'/Sky/')

        





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



'''






















'''


#  SKY BY ELLIPSE FITTING METHOD: TAKES TOO LONG!



# =============================================================================
# Print check
# =============================================================================
print('# =================================================')
print('               Cluster: ' +sim+'-'+cluster+'-'+orientation+'           ')
print('# =================================================')
print('               Fix center: ' +str(fix_center)+'               ')
print('# =================================================')



# =============================================================================
# Elliptical Fitting
# =============================================================================
# Instance our class
elipses = EllipticalProfile(data_masked, Xc, Yc, sma0, ell, theta, fix_center=fix_center)
                           
# Fix the center of the galaxy or Find it. If the algorithm is successful the (x, y) coordinates in this EllipseGeometry 
# (i.e. the x0 and y0 attributes) instance will be modified.

if fix_center == False:
    elipses.find_center()

# We perform the elliptical isophote fit with the fit_image method
elipses.fit_image(sma0=sma0, integrmode='median', step = step, sclip=3, nclip=3, fflag=0.1, maxsma=maxsma)
# We fit the outer regions too if we are going to measure the sky far away from the BCG
#elipses.outer_regions(sma_cut=300, step=2.2*step)


#==============================================================================
# Surface brightness values
#==============================================================================
sbp = elipses.surfbrigth(zp, pixelscale, rms_sky=0, A=0, dimming=0, k_corr=0)



# =============================================================================
# Save our class
# =============================================================================
import pickle
pickle.dump(elipses, open(results_path+'Sky'+'/'+"sky_ellipses.p", "wb" ) )


#==============================================================================
# Read fitted data (from ellipse fitting)
#==============================================================================
r = sbp['sma']

ct = sbp['ct']
ct_err = sbp['ct err']

ell = sbp['ell']
ell_err = sbp['ell err']

pa = sbp['pa']
pa_err = sbp['pa err']

xc = sbp['xc']
xc_err = sbp['xc err']
x_median = np.median(xc)

yc = sbp['yc']
yc_err = sbp['yc err']
y_median = np.median(yc)


ell = sbp['ell']
ell_err = sbp['ell err']        
        
pa = sbp['pa']
pa_err = sbp['pa err']

b4 = sbp['b4']
b4_err = sbp['b4 err']

snr = sbp['SNR']

print('X median: ', x_median)
print('Y median: ', y_median)




# =============================================================================
# Figures to check that all is ok   
# =============================================================================
plt.close('all')


#==============================================================================
# Plot apertures
#==============================================================================
plt.figure(1,figsize=(8,8))
plt.title(sim+'-'+cluster+'-'+orientation, size=25)

plt.imshow(data, origin='lower', cmap='viridis', norm=LogNorm(vmin=0,vmax=10))

from matplotlib import colors
cmap = colors.ListedColormap(['k'])

shadow = np.nan_to_num(mask,nan=1)
shadow[shadow==0]=np.nan

plt.imshow(shadow, cmap=cmap,  interpolation='none', origin='lower',alpha=0.75)

elipses.plot(color = 'r')

plt.xlabel("pix",size=20)
plt.ylabel("pix",size=20)

plt.tight_layout()
plt.savefig(results_path+'Sky'+'/sky_Apertures.pdf')




# =============================================================================
# Profile in counts
# =============================================================================
plt.figure(2,figsize=(8,8)) 
plt.axhline(0,ls='dashed',color='r')

plt.errorbar(r, ct, yerr=ct_err, fmt='o',markersize=4,zorder=0)
plt.xlabel('SMA (pix)',size=14)
plt.ylabel('Counts',size=14)
plt.tight_layout()
plt.show()
plt.savefig(results_path+'Sky'+'/sky_Profile.pdf')



# =============================================================================
# Signal to Noise of each ellipse S/N
# =============================================================================
plt.figure(4)
plt.scatter(r, snr)
plt.axhline(1,ls='dashed',color='C1',zorder=1)
plt.xlabel('SMA (pix)')
plt.ylabel('S/N')
plt.ylim(-0.5,100)
plt.tight_layout()
plt.savefig(results_path+'Sky'+'/sky_SNR.pdf', format='pdf')
# =============================================================================




'''














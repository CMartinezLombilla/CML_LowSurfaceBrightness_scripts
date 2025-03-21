import os, glob, sys, warnings, array, re, math, time, copy
import numpy                  as np
import matplotlib.pyplot      as plt
import matplotlib.cm          as cm
import matplotlib.colors      as mcolors
from   astropy.io             import fits, ascii
from   astropy.wcs            import WCS
from   astropy                import units as u
from   astropy.coordinates    import SkyCoord
from astropy.stats            import sigma_clipped_stats
from astropy.convolution      import Gaussian2DKernel
from astropy.convolution      import convolve
from photutils                import EllipticalAperture
from CML                      import Counts2Sb

plt.close('all')
plt.show()


#######

mag_isophot = np.arange(18.5, 26, 0.5)

zoom_minx = 1000
zoom_maxx = 1500
zoom_miny = 1100
zoom_maxy = 1450

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
    
os.system('rm params.py')


# Def figure
fig = plt.figure(1,figsize=(15,5))


# Images
im_names = ['CoreMasked-g.fits', 'CoreMasked-r.fits', 'CoreMasked-i.fits' ]


## For the reference image :
ref = fits.getdata(results_path+'masks/'+im_names[1])
ref_zoom = ref[zoom_miny:zoom_maxy, zoom_minx:zoom_maxx]

# Counts of each isophote
counts_r = 10**(-(mag_isophot- zp - 5.*np.log10(pixscale) + dimming + k_corrs['r'] + extinctions['r'])/2.5)

# Pixels in each isophote region in this list
regions = []

# Images of the isophotes in each band (should be the same)
isoph_ims = []
isoph_ims_ct = []
sb_mean_from_ct = []


for counts in range(len(counts_r)-1):
    ind_isophot = (ref_zoom < counts_r[counts]) & (ref_zoom > counts_r[counts+1])
    regions.append(ind_isophot)
    
# Load images in each band

for index, im_name in enumerate(im_names):
    image = fits.getdata(results_path+'masks/'+im_name)
    image_zoom = image[zoom_miny:zoom_maxy, zoom_minx:zoom_maxx]
    isoph_im = np.zeros_like(image_zoom)
    isoph_im_ct = np.zeros_like(image_zoom)
    

    band = bands[index]
    sb_im = Counts2Sb(zp, pixscale, A=extinctions[band], dimming=dimming, k_corr=k_corrs[band]).do(image_zoom)
    
    # Create an elliptical aperture around the core to mask beyond it
    # ap = EllipticalAperture((np.shape(image_zoom)[1]/2, 210), a=210, b=130, theta= np.deg2rad(-20))
    # ap_mask = ap.to_mask(method='center').to_image(np.shape(image_zoom))
    # image_zoom[ap_mask==0] = np.nan
    # isoph_im_ct[ap_mask==0] = np.nan
    
    for region in regions:
        median_reg = sigma_clipped_stats(sb_im[region], sigma = 1)[1] # No sig. diff mean-median
        median_reg_ct = sigma_clipped_stats(image_zoom[region], sigma = 1)[1]

        isoph_im[region] = median_reg
        isoph_im_ct[region] = median_reg_ct
        # isoph_im_ct[ap_mask==0] = np.nan

    isoph_ims.append(isoph_im)
    isoph_ims_ct.append(isoph_im_ct)
    
    sb_mean_from_ct.append(Counts2Sb(zp, pixscale, A=extinctions[band], dimming=dimming, k_corr=k_corrs[band]).do(isoph_im_ct))
    
    
    # Plot figures
    ax = plt.subplot(1,3,1+index)
    vmin, vmax = np.percentile(isoph_im, (3, 97))
    isoph_im[isoph_im==0]=np.nan
    ax.imshow(isoph_im, cmap='Blues', origin='lower', interpolation='None')
    

    
# # Saco los colores    
# im_color_gr = isoph_ims[0]-isoph_ims[1]
# im_color_ri = isoph_ims[1]-isoph_ims[2]

    
# Saco los colores    
im_color_gr = sb_mean_from_ct[0]-sb_mean_from_ct[1]
im_color_ri = sb_mean_from_ct[1]-sb_mean_from_ct[2]

color_gr = []
color_ri = []

for region in regions:
    color_gr.append(np.nanmin(im_color_gr[region]))
    color_ri.append(np.nanmin(im_color_ri[region]))
    
    
# Plot color images:    
plt.figure(2)    
plt.imshow(im_color_gr, origin='lower', cmap='YlOrBr', interpolation='None')    
plt.title('g-r color map')
plt.colorbar()

plt.figure(3)    
plt.imshow(im_color_ri, origin='lower', cmap='YlOrBr', interpolation='None')    
plt.title('r-i color map')
plt.colorbar()



# for index, im_name in enumerate(im_names):
#     isoph_im           = np.zeros_like(im)
#     counts = 10**(-(mag_isophot- zp - 5.*np.log10(pixscale) + dimming + k_corrs[bands[index]] + extinctions[bands[index]])/2.5)
    
#     image, hdr = fits.getdata(results_path+'masks/'+im_name, header = True)
#     result      = []

#     for pp in np.arange(0, len(mag_isophot)-1):
#     ## Define la isofota en banda X
#         ind_isophot = (im < counts[pp]) & (im > counts[pp+1])

#         # Generate isophote image
#         isoph_im[ind_isophot] = pp + tt_band
        
        
#         mean_counts = sigma_clipped_stats(image[ind_isophot], sigma = 3)[1]
#         result.append(mean_counts)

#     isophot_counts.append(result)
#     fits.writeto(results_path+'Isophotes/'+'isophotes_'+ bands[tt_band]+'.fits', isoph_im, overwrite = True)

#     # =============================================================================
#     # Zoom in the isophote to the core region
#     # =============================================================================
    
#     isoph_zoom = isoph_im[zoom_min:zoom_max,zoom_min:zoom_max]
    

#     # Figure:
#     #kernel= Gaussian2DKernel(3) #X & Y size by default 8*sigma_pix
#     #conv_isoph_zoom = convolve(isoph_zoom, kernel) 
    
    
    
#     X, Y = np.meshgrid(np.arange(isoph_zoom.shape[0]), np.arange(isoph_zoom.shape[1]))
    
#     ax = plt.subplot(1,3,1+tt_band)
#     ax.contourf(X, Y, isoph_zoom, len(mag_isophot), cmap='Blues', origin='lower')
    
#     # plt.xlim(200,400)
#     # plt.ylim(200,400)
# plt.tight_layout()
    
    
    

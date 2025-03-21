import os, glob, sys, warnings, array, re, math, time, copy
import numpy                  as np
import matplotlib.pyplot      as plt
import matplotlib.cm          as cm
from   astropy.io             import fits, ascii
from   astropy.wcs            import WCS
from   astropy                import units as u
from   astropy.coordinates    import SkyCoord
from   astropy.stats          import sigma_clipped_stats

#######

#redshift      = 0.2
#arcsec_to_kpc = 1./cosmo.arcsec_per_kpc_proper(redshift)
pixscale      = 0.168
arcsec_to_kpc = 1.068 #3.4 # ??


mag_isophot   = np.arange(25, 30, 0.5)

# Aqui es solo una pero vamos a hacer que tenemos varias imagenes
image_name = ['PSFcorr-r-SkySubDist.fits']#, 'PSFcorr-r-SkySubDist.fits' ]
mask_name  = ['CoreMask-r.fits']#, 'CoreMask-r.fits']

## For the reference image :

im, hdr = fits.getdata(image_name[0], header = True)
msk     = fits.getdata(mask_name[0])
im[msk != 0] = np.nan

ZP      = 27#float(hdr['ZP'])
counts  = 10**(-(mag_isophot - ZP)/2.5)

### Define the center of mass from the last isophot.

coords_last = np.argwhere((im < counts[-2]) & (im > counts[-1]))
centerofm   = np.average(coords_last, axis = 0)

isophot_counts = []
test           = np.zeros_like(im)

for tt in np.arange(0, len(image_name)):
    image, hdr  = fits.getdata(image_name[tt], header = True)
    mask        = fits.getdata(mask_name[tt])
    image[mask != 0] = np.nan
    result      = []
    distance    = []
    cm          = []
    for pp in np.arange(0, len(mag_isophot)-1):
        ## Define la isofota en banda X
        ind_isophot = (im < counts[pp]) & (im > counts[pp+1])
        ## Find a "distance"
        pixels      = np.argwhere((im < counts[pp]) & (im > counts[pp+1]))
        dist_vec    = [np.sqrt((pp[0] - centerofm[0])**2  + (pp[1] - centerofm[1])**2)  for pp in pixels]
        distance.append(np.average(dist_vec))
        #Test
        test[ind_isophot] = pp + tt
        mean_counts = sigma_clipped_stats(image[ind_isophot], sigma = 3)[0]
        result.append(mean_counts)

    distance = np.array(distance) * pixscale * arcsec_to_kpc
    isophot_counts.append(result)

    fits.writeto('test_'+ str(tt) + '.fits', test, overwrite = True)


sys.exit()

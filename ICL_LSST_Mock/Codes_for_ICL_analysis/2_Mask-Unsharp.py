#*****************************************************************************
#                       COLD MASK IMAGE PROGRAM   [Rix+2004]
#-----------------------------------------------------------------------------
#                             Cristina Martinez 08/20
#*****************************************************************************

# Versions
# %pip install numpy== 1.18.5
# %pip install matplotlib==3.2.2
# %pip install astropy==4.0
# %pip install photutils==0.7.2

#!pip show numpy
#!pip show matplotlib
#!pip show astropy
#!pip show photutils


# Basics
import numpy as np

# Matplotlib
from astropy.io import fits
import matplotlib.pyplot as plt
plt.close('all')
plt.ion()

# Astropy
from photutils import detect_threshold
from astropy.convolution import Gaussian2DKernel
from photutils import detect_sources
from photutils import deblend_sources
from astropy.convolution import convolve

import os
import glob

from CML import DeleteMaskRegion
from CML import Patchs
from CML import save_fits_image



# =============================================================================
# Initial params
# =============================================================================
#sims = ['Horizon_AGN', 'Hydrangea', 'Magneticum', 'TNG_100'] 
sim = 'Hydrangea' # 'Hydrangea'
clusters = ['299', '338', '353', '727']  
orientations = ['xy','xz', 'yz'] #  
band = 'r'

npix_bin = 1 #number of pixels for the binning
factor = 1e14 # scale factor for multiplying the data values
factor_str = format(factor, '.0e')


# =============================================================================
# Start counter time
# =============================================================================
import time
start = time.time()


# =============================================================================
# Loop for the choosen clusters and orientations
# =============================================================================

for cluster in clusters:    
    for orientation in orientations:
        
        # =============================================================================
        # Paths 
        # =============================================================================
        data_dir = '/Users/z3530031/unsw/Work/Projects/LSST/Chall3-ICL/data/'+sim+'/'
        res_dir = '/Users/z3530031/unsw/Work/Projects/LSST/Chall3-ICL/results/'+sim+'/'+cluster+'_'+orientation+'/'
        
        #image_file = sim+'_'+cluster+'_'+orientation+'_'+band+'_2Mpc'+'_bin'+str(npix_bin)+'_x'+factor_str+'.fits'
        #image_files = glob.glob(image_dir+'*.fits') # Load a list 
        
        
        # =============================================================================
        # Create directory to save results
        # =============================================================================
        try:
            # Create target Directory
            os.mkdir(res_dir+'Masks')
        except:
            pass
        
        
        
        #================================================================================
        # COLD Source detection & deblending -- bright and extended sources
        #================================================================================
        
        # Gaussian kernel for convolve the data for a better source detection
        sigma_pix = 3
        kernel = Gaussian2DKernel(sigma_pix) #X & Y size by default 8*sigma_pix
        kernel.normalize()
        
        
        
        # =============================================================================
        # Unsharp data
        # =============================================================================
        
        # Load images
        im_files = glob.glob(res_dir+'Images/'+'*bin2_x1e+12.fits')
        im, hdr = fits.getdata(im_files[0], header = True) # BUCLE aqui?
        
        # Unsharped image
        # Convolve the image with a gaussian kernel to create an unsharp image
        sigma_pix = 5
        kernel= Gaussian2DKernel(sigma_pix) #X & Y size by default 8*sigma_pix
        conv_im = convolve(im, kernel)
        
        # Subtract data - conv_data
        unsharp = im - conv_im
        
        
        # =============================================================================
        # Multiprocessing of the patchs and source detection/deblending
        # =============================================================================
        
        # Define smaller section over the image as it's is too big
        patchs = Patchs(unsharp, 1715)
        
        
        im_chunks = []
        
        for patch in patchs:
            im_chunks.append(unsharp[patch])
            
            
        
        def patch_mask(im_chunk):
        
            threshold = detect_threshold(im_chunk, nsigma=1.08, sigclip_sigma=2, mask_value = np.nan)
            mask_chunk = np.zeros_like(im_chunk)
            try:
                # Segmentation - detecting sources
                segm = detect_sources(im_chunk, threshold, npixels=60,filter_kernel=kernel)
                # Deblending sources 
                mask_chunk = deblend_sources(im_chunk, segm, npixels=60, filter_kernel=kernel).data
            except:
                print('No object detected in patch ', patch)
            return mask_chunk
        
        
        
        import multiprocessing as mp
        
        pool = mp.Pool(8) 
        results = pool.map(patch_mask, [chunk for chunk in im_chunks])
        
        pool.close()
        
        
        mask = np.zeros_like(im)
        for ind, patch in enumerate(patchs):
            mask[patch] = results[ind]
        
        # =============================================================================
        # End time counter
        # =============================================================================
        print('It took', time.time()-start, 'seconds.')
        
        
        
        
        # =============================================================================
        # Remove mask we don't need over the BCG
        # =============================================================================
        print('Sim: '+sim+'; Cluster: '+cluster+'; Ori: '+orientation)
        print('Remove masks over the BCG')
        
        bcg_patch = patchs[int(len(patchs)/2)]
        mask[bcg_patch] = DeleteMaskRegion(mask[bcg_patch],im[bcg_patch], vmax=2000)
        
        
        
        # =============================================================================
        # Convolve the mask with a Gaussian kernel to make it larger 
        # =============================================================================
        sigma_pix = 3  
        kernel= Gaussian2DKernel(sigma_pix)  #X & Y size by default 8*sigma_pix
        
        conv_mask = convolve(mask, kernel)
        
        # =============================================================================
        # Mask the image with the final mask
        # =============================================================================
        im[conv_mask != 0]= np.nan
        
        
        # =============================================================================
        # Save masked image
        # =============================================================================
        save_fits_image(im, hdr, res_dir+'Masks/', 'unsharp_'+os.path.basename(im_files[0]))

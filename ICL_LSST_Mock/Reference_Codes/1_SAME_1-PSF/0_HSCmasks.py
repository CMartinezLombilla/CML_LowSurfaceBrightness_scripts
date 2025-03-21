#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""""
=============================================================================
    Programa que genera una mascara bool a partir de las regiones 
                    DS9 'bright mask' de HSC
=============================================================================
"""


import matplotlib.pyplot as plt
plt.close('all')

from astropy.io import fits
from astropy.wcs import WCS

from unagi import mask


#==============================================================================
# Parametros iniciales
#==============================================================================
GroupID = 400138
images_dir = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/data/images/'
results_dir = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'
masks_dir = results_dir+str(GroupID)+'/masks/intermediate/'

bands = ['g', 'r', 'i']

#==============================================================================
# Iniciamos diccionarios
#==============================================================================
files = {}
masks = {}
edge = {}
crosstalk = {}
bad = {}
detector_mask = {}

#==============================================================================
# Cargamos las imagenes de las mascaras
#==============================================================================          
for band in bands:
    files[band] = fits.open(images_dir+str(GroupID)+'-'+band+'.fits')
    masks[band] = files[band][2].data
    
    mask_object = mask.Mask(masks[band], data_release='pdr2')

    # =============================================================================
    # SENSOR_EDGE mask -- bordes de los CCDs
    # =============================================================================
    edge[band] = mask_object.extract(['SENSOR_EDGE'])[0]

    # =============================================================================
    # CROSSTALK mask -- reflexiones de estrellas brillantes
    # =============================================================================
    crosstalk[band] = mask_object.extract(['CROSSTALK'])[0]

    # =============================================================================
    # BAD mask -- pixels malos 
    # =============================================================================
    bad[band] = mask_object.extract(['BAD'])[0]

    # =============================================================================
    # Combined mask -- SENSOR_EDGE + CROSSTALK + BAD in each band
    # =============================================================================
    detector_mask[band] = edge[band] + crosstalk[band] + bad[band]

    # =============================================================================
    # Guardamos las nuevas mascaras combinadas 
    # =============================================================================
    detectmask_hdu = fits.PrimaryHDU(detector_mask[band])  # we create a PrimaryHDU object to encapsulate the data:
    file_wcs = files[band][2].header
    wcs = WCS(file_wcs)
    detectmask_hdu.header.update(wcs.to_header())
    detectmask_hdu.header.append(('HISTORY', 'Combined SENSOR_EDGE + CROSSTALK + BAD Mask', 'C. Martinez-Lombilla'), end=True)
    detectmask_hdu.writeto(masks_dir+'detectorHSC-'+band+'.fits', overwrite=True) # write to a new file
    
    
    
plt.imshow(detector_mask[band], cmap='gray', alpha=0.3)


#plt.savefig(results_dir+str(GroupID)+'/GroupBrightObjs.pdf',dpi=300)




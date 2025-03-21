import os
import glob
import sys

import numpy as np

import astropy.units as u
from astropy import wcs
from astropy.wcs import WCS
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.utils.data import download_file, clear_download_cache
from astropy.table import Table

import matplotlib.pyplot as plt

#==============================================================================
# Parametros iniciales
#==============================================================================
groupID_list_yes = [400016,
                    400020,
                    400026,
                    400030,
                    400032,
                    400034,
                    400039,
                    400043,
                    400049,
                    400052,
                    400054,
                    400069,
                    400070,
                    400073,
                    400081,
                    400102,
                    400103,
                    400109,
                    400115,
                    400119,
                    400125,
                    400128,
                    400138,
                    400145,
                    400146,
                    400153,
                    400162,
                    400164,
                    400165,
                    400170,
                    400172,
                    400173,
                    400217,
                    400218,
                    400219,
                    400221,
                    400230,
                    400231,
                    400241,
                    400255,
                    400258,
                    400275,
                    400282,
                    400298,
                    400365,
                    400367,
                    400386,
                    400393]  


images_dir = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/data/images/'
tables_dir = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/data/tables/'


# =============================================================================
# Cargamos las IDs y las separaciones de cada miembro del grupo al centro (IterCen)
# =============================================================================

# Load tables:
tableGal = Table.read(tables_dir+'G3CGalv10.fits')
tableGroup = Table.read(tables_dir+'G3CFoFGroupv10.fits')

# empty lists to save data

groups_params = []

# Get info ra, dec, & separation for all the groups

for GroupID in groupID_list_yes:
    
    catid = tableGal['CATAID'][tableGal['GroupID']==GroupID]
    print(catid)
    
    separations = {}
    for i in catid:
        separations[str(i)] = tableGal['SepIterCen'][tableGal['CATAID']== i][0]
        
    # el tamano de la imagen sera de al menos la separacion maxima al centro del grupo + 20 arcsec   
    cutsize = (max(separations.values()) + 20)* u.arcsec
    
    print(cutsize)
    
    
    
    # Cogemos las coords RA y Dec del grupo de las tablas de GAMA (IterCen)
    
    ra = tableGroup['IterCenRA'][tableGroup['GroupID']==GroupID]
    dec = tableGroup['IterCenDec'][tableGroup['GroupID']==GroupID]
    
    coord_1 = SkyCoord(ra[0], dec[0], frame='icrs', unit='deg')
    
    
    # Descargo y guardo las imagenes
    print('file name: '+str(GroupID)+'-g')

    groups_params.append([coord_1.ra.value, coord_1.dec.value, cutsize.value])

'''

#==============================================================================
# Cargamos unagi para descargar las imagenes
#==============================================================================

from unagi import hsc

# Setup HSC-SSP online data archive
# First, you need to setup a HSC-SSP rerun
# Here we use the ultra-deep field from the PDR2

#Login
#User: lombilla
#Pass: OpeU6cHJm//GrJBd7xpnEXHWET7jwxDTEz6LEYw+


pdr = hsc.Hsc(dr='pdr2', rerun='pdr2_dud')

#pdr = hsc.Hsc(dr='pdr3', rerun='pdr3_dud')



# Design the cutout region
# Using the central coordinate and the size of the cutout region:
#   Central coordinate should be a SkyCoord object in astropy
#   Size of the image should be a Quantity object in astropy. It should have a unit, either angular or physical one is fine
#   For angular size: u.arcsec is default.
#   For physical size: u.kpc is default, and you need to provide the redshift of the object.




urlImage = pdr.form_cutout_url(coord_1, w_half=cutsize*1.8, h_half=cutsize*1.8, filt='HSC-G', mask='on', variance='on')

'''

'''

import webbrowser


webbrowser.open(urlImage) 
input('Press enter when the download is finished --')    

urlImage = pdr2.form_cutout_url(coord_1, w_half=cutsize*1.8, h_half=cutsize*1.8, filt='HSC-R', mask='on', variance='on')
webbrowser.open(urlImage) 
input('Press enter when the download is finished --')    

urlImage = pdr2.form_cutout_url(coord_1, w_half=cutsize*1.8, h_half=cutsize*1.8, filt='HSC-I', mask='on', variance='on')
webbrowser.open(urlImage) 
input('Press enter when the download is finished --') 


#cutImage = hsc_cutout(coord_1, cutout_size=cutsize*1.2, filters='i', archive=pdr2, use_saved=False, output_dir=output_dir, verbose=True, save_output=True)




#==============================================================================
# Combinamos las imagenes gri
#==============================================================================
          
fileG = fits.open(images_dir+str(GroupID)+'-g.fits')
fileR = fits.open(images_dir+str(GroupID)+'-r.fits')
fileI = fits.open(images_dir+str(GroupID)+'-i.fits')

# Data from images
imageG = fileG[1].data
imageR = fileR[1].data
imageI = fileI[1].data

# Header for maintin wcs
wcsG = WCS(fileG[1].header)
wcsR = WCS(fileR[1].header)
wcsI = WCS(fileI[1].header)


#  We sum the 3 images
gri_data = imageG + imageR + imageI


# Save the fits files with wcs
# Create a fits file to save the image 
image_hdu = fits.PrimaryHDU(gri_data)  # we create a PrimaryHDU object to encapsulate the data:
image_hdu.header.update(wcsG.to_header())
image_hdu.header.append(('HISTORY', 'Combined gri', 'C. Martinez-Lombilla - Sum g, r, i bands'), end=True)
image_hdu.writeto(images_dir+str(GroupID)+'-gri.fits', overwrite=True) 

'''
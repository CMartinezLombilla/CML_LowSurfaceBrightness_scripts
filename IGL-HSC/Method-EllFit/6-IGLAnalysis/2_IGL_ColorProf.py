import sys
import os
sys.path.append('/Users/felipe/Dropbox/CML/')
sys.path.append('/Users/z3530031/Dropbox/CML/')

import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
plt.ioff()

from astropy.table import Table


# =============================================================================
# =============================================================================
groupID = '400138'
members = ['IGL']
bands = ['g', 'r', 'i']



# =============================================================================
# Paths
# =============================================================================
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'



# =============================================================================
# Creamos directorios para guardar los distintos productos del programa
# =============================================================================
try:
    # Create target Directory
    os.mkdir(results_path+'IGL/color')
except:
    pass

try:
    # Create target Directory
    os.mkdir(results_path+'IGL/color/profiles')
except:
    pass

# =============================================================================
# Diccionarios para cada uno de los colores
# =============================================================================

colors = ['g-i', 'g-r','r-i']

files = {'g-i':[results_path+'IGL/profiles/IGL_sbp_g.fits', results_path+'IGL/profiles/IGL_sbp_i.fits'],
         'g-r':[results_path+'IGL/profiles/IGL_sbp_g.fits', results_path+'IGL/profiles/IGL_sbp_r.fits'], 
         'r-i':[results_path+'IGL/profiles/IGL_sbp_r.fits', results_path+'IGL/profiles/IGL_sbp_i.fits']}


for color in colors:    
    # =============================================================================
    # Leemos los perfiles que necesitamos: a - b
    # =============================================================================
    
    band_a = files[color][0]
    band_b = files[color][1]

    r = Table.read(band_a)['r'] # La posicion radial e la misma para todos los perfiles
    
    # =============================================================================
    # perfil a
    # =============================================================================
    sb_a  = Table.read(band_a)['sb']
    
    ct_a  = Table.read(band_a)['ct']
    err_a  = Table.read(band_a)['ct err']

    # =============================================================================
    # perfil b
    # =============================================================================
    sb_b  = Table.read(band_b)['sb']
    
    ct_b  = Table.read(band_b)['ct']
    err_b = Table.read(band_b)['ct err']

    
    # =============================================================================
    # color
    # =============================================================================    
    c = sb_a - sb_b
    
    err = -2.5 *np.log10(np.e) * np.sqrt(np.square(err_a/ct_a) + np.square(err_b/ct_b))

    
    file_name = 'IGL_prof_'+color+'.fits'
    table = Table([r,c,err],names=('r','color','err'), meta={'Color Profiles IGL': 'table info'})

    table.write(results_path+'IGL/color/profiles/'+file_name, overwrite=True)
    
    
    # =============================================================================
    # color
    # =============================================================================       
    plt.figure(1)
    plt.errorbar(r, c, yerr=err, marker='.',mew=0, ms=4,linestyle='-',lw=1,elinewidth = 0.5, label = color)
    plt.legend(loc=4, numpoints=1,fontsize=9, framealpha=1)
    plt.xlabel("$R$ [pix]")
    plt.ylabel('Color')
    plt.ylim(-0.1, 1.2)
    plt.xlim(0, 370)


plt.savefig(results_path+'IGL/color/profiles/IGL_ColorProf.pdf')



# -----------------------------------------------------------------------------

# with the galaxies in the image

files = {'g-i':[results_path+'IGL/profiles/group_sbp_g.fits', results_path+'IGL/profiles/group_sbp_i.fits'],
         'g-r':[results_path+'IGL/profiles/group_sbp_g.fits', results_path+'IGL/profiles/group_sbp_r.fits'], 
         'r-i':[results_path+'IGL/profiles/group_sbp_r.fits', results_path+'IGL/profiles/group_sbp_i.fits']}

for color in colors:    
    # =============================================================================
    # Leemos los perfiles que necesitamos: a - b
    # =============================================================================
    
    band_a = files[color][0]
    band_b = files[color][1]

    r = Table.read(band_a)['r'] # La posicion radial e la misma para todos los perfiles
    
    # =============================================================================
    # perfil a
    # =============================================================================
    sb_a  = Table.read(band_a)['sb']
    
    ct_a  = Table.read(band_a)['ct']
    err_a  = Table.read(band_a)['ct err']

    # =============================================================================
    # perfil b
    # =============================================================================
    sb_b  = Table.read(band_b)['sb']
    
    ct_b  = Table.read(band_b)['ct']
    err_b = Table.read(band_b)['ct err']

    # =============================================================================
    # color
    # =============================================================================    
    c = sb_a - sb_b
    
    err = -2.5 *np.log10(np.e) * np.sqrt(np.square(err_a/ct_a) + np.square(err_b/ct_b))

    
    file_name = 'group_prof_'+color+'.fits'
    table = Table([r,c,err],names=('r','color','err'), meta={'Color Profiles IGL+Galaxies': 'table info'})

    table.write(results_path+'IGL/color/profiles/'+file_name, overwrite=True)
    
    
    # =============================================================================
    # color
    # =============================================================================       
    plt.figure(2)   
    plt.errorbar(r, c, yerr=err, marker='.',mew=0, ms=4,linestyle='-',lw=1,elinewidth = 0.5, label = color)
    plt.legend(loc=4, numpoints=1,fontsize=9, framealpha=1)
    plt.xlabel("$R$ [pix]")
    plt.ylabel('Color')
    plt.ylim(-0.1, 1.2)
    plt.xlim(0, 370)


plt.savefig(results_path+'IGL/color/profiles/group_ColorProf.pdf')



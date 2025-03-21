import sys
import os

import numpy as np
import matplotlib.pyplot as plt
plt.ion(),plt.show()
plt.close('all')

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 18

from astropy.table import Table


# Division entera en "base"
def baseround(x, base=10):
    return base * round(x/base)



# =============================================================================
# Parametros iniciales del grupo/imagen/IGL/aperturas
# =============================================================================
groupID = '400138'
members = ['IGL']
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
bcgID = params.bcgID
z = params.z

# Parametros del instrumento
zp = params.zp
pixscale = params.pixelscale


os.system('rm params.py')


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
plot_col = ['C0', 'C1', 'C2']

files = {'g-i':[results_path+'IGL/profiles/dist_IGL_sbp_g.fits', results_path+'IGL/profiles/dist_IGL_sbp_i.fits'],
         'g-r':[results_path+'IGL/profiles/dist_IGL_sbp_g.fits', results_path+'IGL/profiles/dist_IGL_sbp_r.fits'], 
         'r-i':[results_path+'IGL/profiles/dist_IGL_sbp_r.fits', results_path+'IGL/profiles/dist_IGL_sbp_i.fits']}


fig, ax = plt.subplots(figsize=(11,7.5))

for ind, color in enumerate(colors):    
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

    
    file_name = 'dist_IGL_prof_'+color+'.fits'
    table = Table([r,c,err],names=('r','color','err'), meta={'Color Profiles IGL': 'table info'})

    table.write(results_path+'IGL/color/profiles/'+file_name, overwrite=True)
    
    
    # =============================================================================
    # color
    # =============================================================================       
    
    ax.errorbar(r, c, marker='.', color=plot_col[ind],mew=0,ms=4,linestyle='--',lw=1,elinewidth = 0.5)
    ax.fill_between(r, c+err, c-err, color=plot_col[ind], alpha=0.4)
   



#plt.savefig(results_path+'IGL/color/profiles/dist_IGL_ColorProf.pdf')



# -----------------------------------------------------------------------------

# with the galaxies in the image

files = {'g-i':[results_path+'IGL/profiles/dist_group_sbp_g.fits', results_path+'IGL/profiles/dist_group_sbp_i.fits'],
         'g-r':[results_path+'IGL/profiles/dist_group_sbp_g.fits', results_path+'IGL/profiles/dist_group_sbp_r.fits'], 
         'r-i':[results_path+'IGL/profiles/dist_group_sbp_r.fits', results_path+'IGL/profiles/dist_group_sbp_i.fits']}

for ind, color in enumerate(colors):    
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

    
    file_name = 'dist_group_prof_'+color+'.fits'
    table = Table([r,c,err],names=('r','color','err'), meta={'ColorProf IGLandGalx': 'table info'})

    table.write(results_path+'IGL/color/profiles/'+file_name, overwrite=True)
    
    
    # =============================================================================
    # Plot color profiles [+ ticks in kpc]
    # =============================================================================       
    from CML import Pix2Kpc
    from CML import Kpc2Pix

    max_x_pix = 300
    t_l = baseround(Pix2Kpc(max_x_pix, scale=pixscale, z=z))
    delta_t = int(t_l/5)

    
    ticks_labels = np.arange(0,t_l+delta_t,delta_t)
    ticks = Kpc2Pix(ticks_labels, scale=pixscale, z=z)

    
    ax.errorbar(r, c, marker='.', color=plot_col[ind],mew=0, ms=4,linestyle='-',lw=1,elinewidth = 0.5, label = color)
    ax.fill_between(r, c+err, c-err, color=plot_col[ind], alpha=0.6)
    ax.legend(loc='upper right', numpoints=1, framealpha=1)
    ax.set_xlabel("SMA (arcsec)")
    ax.set_ylabel('Color')
    ax.set_ylim(-0.01, 1.3)
    ax.set_xlim(0, max_x_pix)

    ax.yaxis.set_ticks_position('both') 
    ax.xaxis.set_ticks_position('both') 
    
    ax_twin=ax.twiny()
    ax_twin.tick_params(axis="both",direction="in")
    ax_twin.tick_params(axis='x')
    ax_twin.set_xlim(0,max_x_pix)
    ax_twin.set_xticks(ticks)
    ax_twin.set_xticklabels(ticks_labels)
    ax_twin.set_xlabel("SMA (kpc)")

    plt.xlim(0,max_x_pix)





ax2 = ax.twinx()
styles = ['-', '--']
struct = ['Group', 'IGL']

for ss, sty in enumerate(styles):
    ax2.plot(np.NaN, np.NaN, ls=styles[ss],
             label=struct[ss], c='black')
ax2.get_yaxis().set_visible(False)

ax2.legend(loc='upper center')


plt.savefig(results_path+'IGL/color/profiles/dist_group-and-IGL_ColorProf.pdf')



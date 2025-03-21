#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 15:14:50 2021

@author: C. Martinez-Lombilla
"""

import sys
import os
sys.path.append('/Users/felipe/Dropbox/CML/')
sys.path.append('/Users/z3530031/Dropbox/CML/')

import numpy as np
import matplotlib.pyplot as plt
plt.ion(),plt.show()
plt.close('all')


plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 18

plt.rc('text', usetex=False)

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
paper_path = '/Users/z3530031/Dropbox/Aplicaciones/Overleaf/IGL-MNRAS/figures/'


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
sb_lim = params.sb_limit


# A y Correccion K
k_corrs = params.k_corrs
extinctions = params.extinctions
dimming = 10*np.log10(1+z)

os.system('rm params.py')



# =============================================================================
# Diccionarios para cada uno de los ficheros con los perfiles
# =============================================================================


colors = ['C2', 'Crimson', 'C7']
linestyle = {'IGL':'--', 'group': '-'}

fig, ax = plt.subplots(figsize=(9, 6))


for structure in ['IGL', 'group']:

    for index, band in enumerate(bands): 
            
        files = {'IGL':results_path+'IGL/profiles/dist_IGL_sbp_'+band+'.fits',
                 'group': results_path+'IGL/profiles/dist_group_sbp_'+band+'.fits'}
        
        # =============================================================================
        # Leemos los perfiles que necesitamos
        # =============================================================================
    
        prof = Table.read(files[structure])
        r = prof['r']
        sb  = prof['sb']    
        err_up  = prof['err up'] 
        err_down  = prof['err down'] 

      
        # =============================================================================
        # PLOT
        # =============================================================================       
        max_x_pix = 500
        max_x_arcsec = 85
        
        # =============================================================================
        # Ticks en Kpc
        # =============================================================================
        from CML import Pix2Kpc
        from CML import Kpc2Pix
        t_l = baseround(Pix2Kpc(max_x_pix, scale=pixscale, z=z))
        delta_t = int(t_l/5)
    
        
        ticks_labels = np.arange(0,t_l+delta_t,delta_t)
        ticks = Kpc2Pix(ticks_labels, scale=pixscale, z=z)
    
    
        # ax.errorbar(r, sb_IGL, yerr=[err_IGL_down,err_IGL_up],color=colors[index], marker='o',mew=0, ms=4,linestyle=':',lw=2,elinewidth = 0.5, label = 'IGL '+band)
        #plt.errorbar(r, sb_gal, yerr=[err_gal_down,err_gal_up],color=colors[index], marker='s',mew=0, ms=4,linestyle='--',lw=1,elinewidth = 0.5, label = 'IGL+Gal '+band)
    
        #plt.plot(r, sb_gal, color=colors[index],linestyle='-',lw=2, label = 'IGL+Gal '+band)
        ax.plot(r*pixscale, sb, color=colors[index], linestyle=linestyle[structure], label=band)
        ax.fill_between(r*pixscale, sb+err_up, sb-err_down, color=colors[index], alpha=0.4)
        
        # Add horixontal lines with thw sb limits corrected by sb-dimming & k-corr
        sb_lim_corr = sb_lim[band] - dimming - k_corrs[band] 
        ax.hlines(sb_lim_corr, 0, r[-1]*pixscale, color=colors[index], linestyle=(0, (1, 3)), zorder=10, alpha=0.8)
        #ax.grid('dashed')
        

        ax.set_ylim(30,18.5)
        ax.set_xlim(0,max_x_arcsec)
        ax.set_xlabel('SMA (arcsec)')
        ax.set_ylabel('$\mu$ (mag/arcsec$^{2}$)')
        ax.tick_params(left=True) 
        ax.yaxis.get_ticklocs(minor=True)
        ax.minorticks_on()
        ax.tick_params(left=True, which='both', direction="out")
     
        
# =============================================================================
#  Ejes en Kpc
# =============================================================================
ax_twin=ax.twiny()
ax_twin.yaxis.get_ticklocs(minor=True)
ax_twin.minorticks_on()
ax_twin.tick_params(axis='x', which='both', direction="out")
ax_twin.set_xlim(0,max_x_pix)
ax_twin.set_xticks(ticks)
ax_twin.set_xticklabels(ticks_labels)
ax_twin.set_xlabel('SMA (kpc)')
plt.xlim(0,max_x_pix)

ay_twin=ax.twinx()
ay_twin.yaxis.get_ticklocs(minor=True)
ay_twin.minorticks_on()
ay_twin.tick_params(axis="both",  which='both', direction="in")
ay_twin.yaxis.set_ticklabels([])
ay_twin.set_ylim(30,18.5)
        
# HSC FWHM R~0.7-8 arcsec r-band 
ax.axvspan(-0.01, 0.8, alpha=0.2, color='C0', zorder=9)  
ax.text(2, 28, 'HSC FWHM r-band', rotation='vertical', size=12, color='C0') 
   
handles, labels = ax.get_legend_handles_labels()

ax.legend(handles[3:], labels[3:], numpoints=1, framealpha=1, loc='upper right')

ax2 = ax.twinx()

styles = ['-', '--']
struct = ['Group', 'IGL']

for ss, sty in enumerate(styles):
    ax2.plot(np.NaN, np.NaN, ls=styles[ss],
             label=struct[ss], c='black')
ax2.get_yaxis().set_visible(False)

ax2.legend(loc='upper center')

plt.savefig(paper_path+'/IGL_sbp.pdf')



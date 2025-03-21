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
# Diccionarios para cada uno de los ficheros con los perfiles
# =============================================================================


colors = ['C8', 'C6', 'C7']
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
        max_x = 500
        
        # =============================================================================
        # Ticks en Kpc
        # =============================================================================
        from CML import Pix2Kpc
        from CML import Kpc2Pix
        t_l = baseround(Pix2Kpc(max_x, scale=pixscale, z=z))
        delta_t = int(t_l/5)
    
        
        ticks_labels = np.arange(0,t_l+delta_t,delta_t)
        ticks = Kpc2Pix(ticks_labels, scale=pixscale, z=z)
    
    
        # ax.errorbar(r, sb_IGL, yerr=[err_IGL_down,err_IGL_up],color=colors[index], marker='o',mew=0, ms=4,linestyle=':',lw=2,elinewidth = 0.5, label = 'IGL '+band)
        #plt.errorbar(r, sb_gal, yerr=[err_gal_down,err_gal_up],color=colors[index], marker='s',mew=0, ms=4,linestyle='--',lw=1,elinewidth = 0.5, label = 'IGL+Gal '+band)
    
        #plt.plot(r, sb_gal, color=colors[index],linestyle='-',lw=2, label = 'IGL+Gal '+band)
        ax.plot(r, sb, color=colors[index], linestyle=linestyle[structure], label=band)
        ax.fill_between(r, sb+err_up, sb-err_down, color=colors[index], alpha=0.4)
    
        #ax.grid('dashed')
        
        ax.set_ylim(30,18.5)
        ax.set_xlim(0,max_x)
        
        ax.yaxis.set_ticks_position('both') 
        ax.xaxis.set_ticks_position('both') 
        ax.set_xlabel('SMA (pix)')
        
        ax.set_ylabel('$\mu$ (mag/arcsec$^{2}$)')
            
     
        
        # =============================================================================
        #  Ejes en Kpc
        # =============================================================================
        ax_twin=ax.twiny()
        ax_twin.tick_params(axis="both",direction="in")
        ax_twin.tick_params(axis='x')
        ax_twin.set_xlim(0,max_x)
        ax_twin.set_xticks(ticks)
        ax_twin.set_xticklabels(ticks_labels)
        ax_twin.set_xlabel("SMA (kpc)")
        plt.xlim(0,max_x)

    
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

plt.savefig(results_path+'IGL/profiles/dist_Group-and-IGL_prof-AllBands.pdf')



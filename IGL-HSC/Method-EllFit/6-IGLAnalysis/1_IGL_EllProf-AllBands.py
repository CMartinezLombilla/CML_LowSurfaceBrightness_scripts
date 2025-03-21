#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 15:14:50 2021

@author: C. Martinez-Lombilla

Choose the structure to get the profiles: "Group" (this also includes a profile of the IGL) or "IGL only" 

"""

import sys
import os
sys.path.append('/Users/felipe/Dropbox/CML/')
sys.path.append('/Users/z3530031/Dropbox/CML/')

import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14

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


files = {'g':[results_path+'IGL/profiles/IGL_sbp_g.fits', results_path+'IGL/profiles/group_sbp_g.fits'],
         'r':[results_path+'IGL/profiles/IGL_sbp_r.fits', results_path+'IGL/profiles/group_sbp_r.fits'], 
         'i':[results_path+'IGL/profiles/IGL_sbp_i.fits', results_path+'IGL/profiles/group_sbp_i.fits']}

colors = ['C6', 'C8', 'C7']




for structure in ['group', 'IGL']:
    fig,ax = plt.subplots(figsize=(9, 5))
    
    for index, band in enumerate(bands):
        
        # =============================================================================
        # Leemos los perfiles que necesitamos: a - b
        # =============================================================================
        
        profiles = {'IGL': results_path+'IGL/profiles/IGL_sbp_'+band+'.fits', 
                'group': results_path+'IGL/profiles/group_sbp_'+band+'.fits'}
        
    
        # =============================================================================
        # Cargo los perfiles
        # =============================================================================
        r = Table.read(profiles['IGL'])['r'] # La posicion radial es la misma para todos los perfiles
    
        sb_IGL  = Table.read(profiles['IGL'])['sb']    
        err_IGL_up  = Table.read(profiles['IGL'])['err up'] 
        err_IGL_down  = Table.read(profiles['IGL'])['err down'] 
    
        sb_gal  = Table.read(profiles['group'])['sb']
        err_gal_up  = Table.read(profiles['group'])['err up'] 
        err_gal_down  = Table.read(profiles['group'])['err down'] 
    
      
        # =============================================================================
        # PLOT
        # =============================================================================       
        max_x = 400
        
        # =============================================================================
        # Ticks en Kpc
        # =============================================================================
        from CML import Pix2Kpc
        from CML import Kpc2Pix
        t_l = baseround(Pix2Kpc(max_x,scale=pixscale,z=z))
        delta_t = int(t_l/5)
    
        
        ticks_labels = np.arange(0,t_l+delta_t,delta_t)
        ticks = Kpc2Pix(ticks_labels, scale=pixscale,z=z)
        
        if structure == 'group':
            # Group + IGL profiles
            ax.plot(r, sb_gal, color=colors[index],linestyle=':',lw=2, label = 'Group '+band)
            ax.fill_between(r, sb_gal+err_gal_up, sb_gal-err_gal_down, color=colors[index], alpha=0.4)
        
        # IGL only 
        ax.plot(r, sb_IGL, color=colors[index],lw=1, label = 'IGL '+band)
        ax.fill_between(r, sb_IGL+err_IGL_up, sb_IGL-err_IGL_down, color=colors[index], alpha=0.4)
        
        # Common set ups
        ax.grid('dashed')
        ax.legend(numpoints=1,fontsize=12, framealpha=1)
        ax.set_ylim(31,21.5)
        ax.set_xlim(0,max_x)
        
        ax.yaxis.set_ticks_position('both') 
        ax.xaxis.set_ticks_position('both') 
        ax.set_xlabel('SMA (pix)',size=14)
        
        ax.set_ylabel('$\mu_{'+band+'}$ (mag/arsec$^{2}$)',size=14)
            
     
        
        # =============================================================================
        #  Ejes en Kpc
        # =============================================================================
        ax_twin=ax.twiny()
        ax_twin.tick_params(axis="both",direction="in")
        ax_twin.tick_params(axis='x')
        ax_twin.set_xlim(0,max_x)
        ax_twin.set_xticks(ticks)
        ax_twin.set_xticklabels(ticks_labels)
        ax_twin.set_xlabel("SMA (kpc)",size=14)
        plt.xlim(0,max_x)
    
    
    plt.savefig(results_path+'IGL/profiles/'+structure+'_prof-AllBands.pdf')
    



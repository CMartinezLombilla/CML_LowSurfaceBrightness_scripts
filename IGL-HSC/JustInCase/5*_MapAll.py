import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from adjustText import adjust_text # conda install -c conda-forge adjusttext 


plt.close('all')
plt.ioff()

plt.rcParams["font.family"] = "serif"
plt.rc('text', usetex=False)

from astropy.io import fits

# =============================================================================
# importamos CML
# =============================================================================

import sys
sys.path.append('/Users/felipe/Dropbox/CML/') # path del CML Feli
sys.path.append('/Users/z3530031/Dropbox/CML/') # path del CML Cris

sys.path.append('/Users/felipe/Dropbox/CML/WHT/') # path del fichero de parametros Feli
sys.path.append('/Users/z3530031/Dropbox/CML/WHT/') # path del fichero de parametros Cris

from CML import LogNorm

#==============================================================================
# Parametros iniciales
#==============================================================================
GroupID = '400282'

# Imagen a utilizar (guardade en  groupID+'/0_data/'):
    
im_file ="gri-cut.fits"
# =============================================================================
# Fichero de parametros
# =============================================================================
# Copiamos el fichero de parametros para asegurarnos de que es el mismo

cmd = 'cp '+GroupID+'/params.py .'

os.system(cmd)
import params

# Parametros del grupo
bcgID = params.bcgID
z = params.z

# Parametros del instrumento
zp = params.zp
pixelscale = params.pixelscale



# Parametros de la galaxia
Xc = params.Xc
Yc = params.Yc

os.system('rm params.py')


members = [*Xc][1:]


# =============================================================================
# Leemos la mascara, los datos y los enmascaramo 
# =============================================================================

mask = fits.getdata(GroupID+'/0_Masks/GroupMask.fits')
data = fits.getdata(GroupID+'/0_data/'+im_file)




# =============================================================================
# Figuras de comprobacion        
# =============================================================================
plt.close('all')
#==============================================================================
# Comprobacion: plot de las aperturas
#==============================================================================


plt.figure(1,figsize=(8,8))
plt.title('Group ID: '+GroupID,size=25)

# =============================================================================
# Datos sin enmascarar
# =============================================================================
plt.imshow(data, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=10))


# =============================================================================
# Sombreamos las zonas enmascaradas
# =============================================================================
from matplotlib import colors
cmap = colors.ListedColormap(['k'])

shadow = np.nan_to_num(mask,nan=1)
shadow[shadow==0]=np.nan

plt.imshow(shadow,cmap=cmap, interpolation='nearest',origin='lower',alpha=0.75)


# =============================================================================
# Nombre de cada uno de los miembros y submiembros (i.e. Gal# y  a,b,c...)
# =============================================================================
# Hacemos unas listas con las etiquetas de miembros y submiembros con sus coordenadas. asignamos a
label = []
x = []
y = []

name= ''
for i in members:
    n = i.split('-')   
    if name != n[0]:
        name = n[0]
        label.append(name)
        x.append(Xc[i])
        y.append(Yc[i])
    
    # Pinta redondeles en todas las posiciones
    if len(n)==1:
        plt.scatter(Xc[i],Yc[i], facecolors='none', s=10,color='w',lw=0.4)

    if len(n)>1:
        label.append(n[1])
        x.append(Xc[i])
        y.append(Yc[i])
        plt.scatter(Xc[i],Yc[i], facecolors='none', s=10,color='w',lw=0.4)

texts = [plt.text(x[i], y[i],label[i], va='center',ha='center',color='w') for i in range(len(label))]
adjust_text(texts)


plt.xticks([],[])
plt.yticks([],[])
plt.tight_layout()

plt.savefig(GroupID+'/Map.pdf',dpi=300)

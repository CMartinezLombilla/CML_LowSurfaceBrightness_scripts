# =============================================================================
# Basics
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
# =============================================================================
# Astropy
# =============================================================================
from astropy.io import fits
from astropy.wcs import WCS
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip

# =============================================================================
# importamos CML
# =============================================================================

import sys
sys.path.append('/Users/felipe/Dropbox/CML/') # path del CML Feli
sys.path.append('/Users/z3530031/Dropbox/CML/') # path del CML Cris

sys.path.append('/Users/felipe/Dropbox/CML/WHT/') # path del fichero de parametros Feli
sys.path.append('/Users/z3530031/Dropbox/CML/WHT/') # path del fichero de parametros Cris

from CML import Centroid
from CML import LogNorm
from CML import LinePhot


#==============================================================================
#Initial parameters
#==============================================================================
groupID = '400138'

member = 'BCG'

# =============================================================================
# Centros aproximados 
# =============================================================================

center_0 ={'1660730': [768, 394],
           '1660646': [816, 385],
           'BGG': [844, 352]}

# radio de busqueda del centroide
r ={'1660730': 6,
    '1660646': 6,
    'BGG': 6}

# Anchura y altuda del nuevo muestreo de la linea que une ambos centros
w = 1
h = 1.5



keys = [*center_0]

# =============================================================================
# PDF para las imagenes de control
# =============================================================================
pdf = PdfPages(groupID+'/0_Masks/'+member+'_Mask.pdf')

# =============================================================================
# Grupos dos a dos
# =============================================================================
n_group = len(keys)-1

groups = []

for i in range(n_group): 
    groups.append([keys[i],keys[i+1]])
    

# =============================================================================
# Leemos la mascara, los datos y los enmascaramo 
# =============================================================================
mask = fits.getdata(groupID+'/0_Masks/GroupMask.fits')    
data = fits.getdata(groupID+'/0_data/gri-cut.fits')

wcs = WCS(fits.getheader(groupID+'/0_data/gri-cut.fits'))

# =============================================================================
#  Ajuste de los centros de todos lo centros
# =============================================================================
xc = {}
yc = {}

for key in keys:
    xc[key],yc[key] = Centroid(data,np.array(center_0[key]),r=r[key])



mask_2 = np.copy(mask)

# =============================================================================
# Loop que recorre todos los grupos
# =============================================================================
for group in groups:
    mask_1 = np.copy(mask_2)


    k1 = group[0]
    k2 = group[1]
    
    x_1 = xc[k1]
    x_2 = xc[k2]
    
    y_1 = yc[k1]
    y_2 = yc[k2]
    
    # =============================================================================
    # Aperturas en la linea que une ambos centros
    # =============================================================================
    # Los limites seran los picos del ajuste
    xmin = np.min([x_1,x_2])
    xmax = np.max([x_1,x_2])
    
    # Intanciamos la clase que nos hace el mustreo a lo largo de la recta que une los centros
    recta = LinePhot(x_1, x_2, y_1, y_2, xmin, xmax, w, h)

    # =============================================================================
    # Ajuste del corte
    # ============================================================================= 
    # Cuentas a lo largo de de la recta que uno los centros
    s, c =recta.phot(data, method='center')
    
    # Polinomio de grado 2 con sigma clip
    p_init =  models.Polynomial1D(degree=2)
    fit = fitting.LinearLSQFitter()
    sigma_fit = fitting.FittingWithOutlierRemoval(fit, sigma_clip, niter=3, sigma=2.0)
    p = sigma_fit(p_init, s, c)[0]
    
    # Vertice de la parabola
    b = p.c1.value
    a = p.c2.value
    
    # EL punto determinado por el vertice de la parabola
    Xv =-b/(2*a)
    Yv = recta.centerLine(Xv)

    # =============================================================================
    # Creamos las dos mascaras
    # =============================================================================

    x_dir = range(mask.shape[1])
    y_dir = range(mask.shape[0])
    
    
    for i in x_dir:
        for j in y_dir:
            if j <= recta.bondary(i,Xv,Yv):
                mask_1[j,i]=np.nan
            else:
                mask_2[j,i]=np.nan
    
    # =============================================================================
    # Guardamos las dos mascaras salvando con la misma cabecera que la original
    # =============================================================================
    
    # Create a fits file to save the first  mask 
    mask_1_hdu = fits.PrimaryHDU(mask_1)  # we create a PrimaryHDU object to encapsulate the data:
    mask_1_hdu.header.update(wcs.to_header())
    mask_1_hdu.writeto(groupID+'/0_Masks/'+member+'-'+k1+'GroupMask.fits', overwrite=True) # write to a new file
    
    # =============================================================================
    # Grafico de control
    # =============================================================================
    plt.close('all')
    fig = plt.figure(1,figsize=(10,5))
    plt.suptitle(k1+'-'+k2+' Mask', size=16)
    # =============================================================================
    #  Marcas de las posiciones ajustadas y las aperturas definidas
    # =============================================================================
    ax1 = fig.add_subplot(121)
    ax1.yaxis.set_ticks_position('both') 
    ax1.xaxis.set_ticks_position('both') 
    ax1.tick_params(direction='in',which='both')  
    
    plt.imshow(data, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=40))

    # Plot de las aperturas definidas
    recta.plot(color='k',ls='dashed')
    
    plt.scatter(x_1,y_1,marker='+',color='r',zorder=2)
    plt.scatter(x_2,y_2,marker='+',color='b',zorder=2)
    
    # distancia entre los puntos ajustados
    d = np.sqrt(np.sum(np.square(np.array([x_1,y_1])-np.array([x_2,y_2]))))
    
    xlim = (Xv-d,Xv+d)
    ylim = (Yv-d,Yv+d)
    
    x =np.arange(Xv-d,Xv+d)
    
    y = recta.bondary(x, Xv,Yv)
    
    # Recta limitrofe entre las mascaras
    plt.plot(x,y,color='r',ls='dashed')
    
    
    plt.xlabel('X [pix]',size=14)
    plt.ylabel('Y [pix]',size=14)
    
    plt.xlim(xlim)
    plt.ylim(ylim)
    
    
    # =============================================================================
    #  Ajust de la posicion
    # =============================================================================
    ax2 = fig.add_subplot(122)
    
    plt.plot(s,c,'o', color='C0')
    
    # Pintamos el ajuste con mas muestreo
    sf =np.linspace(s.min(),s.max(),100)
    plt.plot(sf,p(sf),lw=1.2, color='C2',ls='dashed')
    
    plt.axvline(Xv,color='r',ls='dashed')
    
    plt.xlabel('Proyected coordinates [pix]',size=14)
    plt.ylabel('Counts',size=14)
    
    plt.tight_layout()
    
    plt.show()
    
    pdf.savefig(1)

# Create a fits file to save the  last mask 
mask_2_hdu = fits.PrimaryHDU(mask_2)  # we create a PrimaryHDU object to encapsulate the data:
mask_2_hdu.header.update(wcs.to_header())
mask_2_hdu.writeto(groupID+'/0_Masks/'+member+'-'+k2+'GroupMask.fits', overwrite=True) # write to a new file
    
pdf.close()













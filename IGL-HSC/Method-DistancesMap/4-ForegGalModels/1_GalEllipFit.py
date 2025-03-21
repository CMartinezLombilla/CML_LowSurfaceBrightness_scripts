import numpy as np
import numpy.ma as ma
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.close('all')
plt.ioff()

plt.rcParams["font.family"] = "serif"
plt.rc('text', usetex=False)

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from photutils.aperture import EllipticalAperture

# =============================================================================
# importamos CML
# =============================================================================

import sys
sys.path.append('/Users/felipe/Dropbox/CML/') # path del CML Feli
sys.path.append('/Users/z3530031/Dropbox/CML/') # path del CML Cris

sys.path.append('/Users/felipe/Dropbox/CML/WHT/') # path del fichero de parametros Feli
sys.path.append('/Users/z3530031/Dropbox/CML/WHT/') # path del fichero de parametros Cris

from CML import EllipticalProfile
from CML import LogNorm

# Division entera en "base"
def baseround(x, base=10):
    return base * round(x/base)


# =============================================================================
# Parametros de las galaxias a ajustar
# =============================================================================
groupID = '400138'
bands = ['g', 'r', 'i']

# CENTROS
n_galToSub = 2
Xc_dic = {'g1':1264.603524860944, 'g2':1272.5763834652637} #, 'g3':1213.9558374702103, 'g4':1195.4976329990789}
Yc_dic = {'g1':1233.2791329942963, 'g2':1215.0445739902905} #, 'g3':1241.223156521232, 'g4':1266.00}

# Parametros iniciales del ellipse fitting 
fix_center = True
sma0 = [8, 4]
minsma= [2, 2]
eps = [0.3, 0.1]
pa = [1e-6, 1e-6]
step = [0.2, 0.2] 
maxsma = [30, 15]


# =============================================================================
# Tamano de la region alrededor del core a plotear (radio en pix)
# =============================================================================
zoom_r = 200


# =============================================================================
# Paths
# =============================================================================
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'
masks_path = results_path+'masks/'


# =============================================================================
# Creamos un directorio para guardar los distintos productos del programa
# =============================================================================
try:
    # Create target Directory
    os.mkdir(results_path+'GalaxySubtraction')
except:
    pass

try:
    # Create target Directory
    os.mkdir(results_path+'GalaxySubtraction/profiles')
except:
    pass




for band in bands:

    # =============================================================================
    # Comenzamos cargando los datos y mascaras necesarias
    # =============================================================================
    
    mask_file = masks_path+'intermediate/GroupMask-'+band+'.fits'
    mask = fits.getdata(mask_file)
    ellMask = np.zeros_like(mask)
 
    
    for gal in range(n_galToSub):
 
        gal_number = gal+1
            

        if gal_number==1:
            data_file = results_path+'PSF/PSFcorr-'+band+'-SkySubDist.fits'
        # Second star and onwards
        else: 
            data_file = results_path+'GalaxySubtraction/'+band+'-gal'+str(gal_number-1)+'.fits'    
        
        # =============================================================================
        # Leemos la mascara, los datos y los enmascaramos
        # =============================================================================
        
        data = fits.getdata(data_file)
        header = fits.getheader(data_file)
        
        
        shape = data.shape
        data_masked  = ma.masked_array(data, mask)
        
        
        # =============================================================================
        # Centros
        # =============================================================================
        '''
        # Calculamos Xc e Yc de la posicion de la estrella
        coord = SkyCoord(ra,dec,unit='deg')
        wcs = WCS(header)
        
        Xc,Yc = wcs.world_to_pixel(coord)
        '''
        Xc, Yc = Xc_dic['g'+str(gal_number)], Yc_dic['g'+str(gal_number)]
        
        # =============================================================================
        # Elliptical Fitting
        # =============================================================================
        
        # Instaciamos nuestra clase
        elipses = EllipticalProfile(data, Xc, Yc, sma0[gal], eps=eps[gal], pa=pa[gal], fix_center=fix_center, fix_pa=False, fix_eps=False)
        
        # Si no está fijo el centro lo ajustamos                            
        # Find the center of a galaxy. If the algorithm is successful the (x, y) coordinates in this EllipseGeometry (i.e. the x0 and y0 attributes) instance will be modified.
        
        if fix_center != True:
            elipses.find_center()
        
        
        # We perform the elliptical isophote fit with the fit_image method
        elipses.fit_image(sma0=sma0[gal], minsma=minsma[gal], maxsma=maxsma[gal], step=step[gal], integrmode='median', sclip=3., nclip=3, fflag=0.7)
        
        # =============================================================================
        # Guardamos la clase instanciada -- el ajuste!
        # =============================================================================
        import pickle
        
        pickle.dump(elipses, open(results_path+'GalaxySubtraction/profiles/'+band+'-gal'+str(gal_number)+'ellFit.p', 'wb' ) )

        #==============================================================================
        # leemos los datos ajustados
        #==============================================================================
        
        r = elipses.isolist.sma
        """
        sb = sbp['sb']
        err_up = sbp['err up']
        err_down = sbp['err down']
        
        ct = sbp['ct']
        ct_err = sbp['ct err']
        
        
        ell = sbp['ell']
        ell_err = sbp['ell err']
        
        pa = sbp['pa']
        pa_err = sbp['pa err']
        """
        xc = elipses.isolist.x0
        xc_err = elipses.isolist.x0_err
        x_median = np.median(xc)
        
        yc = elipses.isolist.y0
        yc_err = elipses.isolist.y0_err
        y_median = np.median(yc)
        
        """
        ell = sbp['ell']
        ell_err = sbp['ell err']        
                
        pa = sbp['pa']
        pa_err = sbp['pa err']
        
        b4 = sbp['b4']
        b4_err = sbp['b4 err']
        
        snr = sbp['SNR']
        """
        print('X median: ', x_median)
        print('Y median: ', y_median)
        
        
        # =============================================================================
        # 2D model     
        # =============================================================================
        
        from photutils.isophote import build_ellipse_model
        model_image = build_ellipse_model(data.shape, elipses.isolist[2:-2])
        residual = data - model_image
        
        
        
        # =============================================================================
        # Figuras de comprobacion        
        # =============================================================================
        plt.close('all')
        # =============================================================================
        # PDF para las imagenes de control
        # =============================================================================
        pdf = PdfPages(results_path+'GalaxySubtraction/'+band+'-gal'+str(gal_number)+'_Subtraction.pdf')
                       
        #==============================================================================
        # Comprobacion: plot de las aperturas
        #==============================================================================
        xlim = np.shape(data_masked)[0]/2-zoom_r, np.shape(data_masked)[0]/2+zoom_r
        ylim = np.shape(data_masked)[1]/2-zoom_r, np.shape(data_masked)[1]/2+zoom_r
        
        plt.figure(1,figsize=(8,8))
        plt.title('Group ID: '+groupID,size=25)
        # Datos sin enmascarar
        plt.imshow(data, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=10))
        
        
        # =============================================================================
        # Sombreamos las zonas enmascaradas
        # =============================================================================
        from matplotlib import colors
        cmap = colors.ListedColormap(['k'])
        
        shadow = np.nan_to_num(mask,nan=1)
        shadow[shadow==0]=np.nan
        
        plt.imshow(shadow,cmap=cmap,interpolation='nearest',origin='lower',alpha=0.65)
        
        # Aperturas
        elipses.plot(color = 'r')
        plt.xlim(xlim)
        plt.ylim(ylim)
        
        
        
        # more plots
        plt.figure(4,figsize=(8,8))
        plt.title('Group ID: '+groupID,size=25)
        plt.imshow(residual, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=10))
        plt.xlim(xlim)
        plt.ylim(ylim)
        
        plt.figure(3,figsize=(8,8))
        plt.title('Group ID: '+groupID,size=25)
        plt.imshow(data, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=10))
        plt.xlim(xlim)
        plt.ylim(ylim)
        
        plt.figure(2,figsize=(8,8))
        plt.title('Group ID: '+groupID,size=25)
        plt.imshow(model_image, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=10))
        plt.xlim(xlim)
        plt.ylim(ylim)
        
        plt.show()
        
        for fig in range(1, plt.figure().number): ## will open an empty extra figure :(
            pdf.savefig( fig )
        pdf.close()
        
        
        # =============================================================================
        # Save image with the galaxy subtracted
        # =============================================================================
        
        res = fits.PrimaryHDU(residual)  # we create a PrimaryHDU object to encapsulate the data:
        wcs = WCS(header)
        res.header.update(wcs.to_header())
        res.writeto(results_path+'GalaxySubtraction/'+band+'-gal'+str(gal_number)+'.fits', overwrite=True) # write to a new file
          
        
        # =============================================================================
        # Guardamos las mascaras de cada galaxia cogiendo la penultima ellipse                 
        # =============================================================================
            
        a = elipses.isolist[-2]
        
        ellAp = EllipticalAperture([a.x0, a.y0], a.sma, a.sma*(1 - a.eps), a.pa)
        
        ellMask = ellMask + ellAp.to_mask('center').to_image(shape)
        
        
           
    # =============================================================================
    # Guardamos la mascara con las galaxias sustraidas enmascaradas
    # =============================================================================
    
    mask[ellMask!=0]=np.nan
    
    new_mask = fits.PrimaryHDU(mask)  # we create a PrimaryHDU object to encapsulate the data:
    wcs = WCS(header)
    new_mask.header.update(wcs.to_header())
    new_mask.writeto( masks_path+'GroupMask-'+band+'.fits', overwrite=True) # write to a new file
        
    
    






















"""
# =============================================================================
# Comprobacion: plot de los parametros
# =============================================================================

plt.figure(2, figsize=(8, 8))
plt.suptitle('Group ID: '+GroupID+' - CATAID: '+member,size=25)

plt.subplots_adjust(hspace=0.35, wspace=0.35)

# Ticks y etiquetas para el resto de figuras han de ser los mismo que en la figura 2D pero en la parte positiva

ticks = ticks[[ticks_labels>=0]]
ticks_labels = ticks_labels[ticks_labels>=0]



# =============================================================================
# Ellipticity
# =============================================================================
ax0 = plt.subplot(2, 2, 1)
plt.plot(r, ell,linestyle='solid')
plt.fill_between(r, ell+ell_err, ell-ell_err, alpha=0.4)
plt.xlim(0,maxsma)
#plt.xscale('symlog',subsx = range(2,10),linthresh=1e-1)
#plt.xticks([0,1,10,100],[0,1,10,100])
plt.xlabel('SMA (pix)',size=14)
plt.ylabel('Ellipticity',size=14)

# =============================================================================
#  Ejes en Kpc
# =============================================================================
ax_twin=ax0.twiny()
ax_twin.tick_params(axis="both",direction="in")
ax_twin.tick_params(axis='x')
ax_twin.set_xlim(0,maxsma)
ax_twin.set_xticks(ticks)
ax_twin.set_xticklabels(ticks_labels)
ax_twin.set_xlabel("SMA (kpc)",size=14)
plt.xlim(0,maxsma)
# =============================================================================


# =============================================================================
# P. Angle
# =============================================================================
ax0 = plt.subplot(2, 2, 2)
plt.plot(r, pa/np.pi*180.,linestyle='solid')
plt.fill_between(r, pa/np.pi*180.+pa_err/np.pi*180., pa/np.pi*180.-pa_err/np.pi*180., alpha=0.4)

plt.xlim(0,maxsma)
#plt.xscale('symlog',subsx = range(2,10))
#plt.xticks([0,1,10,100],[0,1,10,100])
plt.xlabel('SMA (pix)',size=14)
plt.ylabel('PA (deg)',size=14)

# =============================================================================
#  Ejes en Kpc
# =============================================================================
ax_twin=ax0.twiny()
ax_twin.tick_params(axis="both",direction="in")
ax_twin.tick_params(axis='x')
ax_twin.set_xlim(0,maxsma)
ax_twin.set_xticks(ticks)
ax_twin.set_xticklabels(ticks_labels)
ax_twin.set_xlabel("SMA (kpc)",size=14)
plt.xlim(0,maxsma)
# =============================================================================

# =============================================================================
# Xc, Yc 
# =============================================================================
ax0 = plt.subplot(2, 2, 3)
plt.errorbar(r, xc, color='C0', fmt='o',markersize=4,zorder=0)
plt.errorbar(r, yc, color='C0', fmt='o',markersize=4,zorder=0)

plt.axhline(x_median,ls='dashed',color='C1',label='median X0: '+str(round(x_median, 2)),zorder=1)
plt.axhline(y_median,ls='dashed',color='C2',label='median Y0: '+str(round(y_median, 2)),zorder=1)

plt.xlim(0,maxsma)
#plt.xscale('symlog',subsx = range(2,10))
#plt.xticks([0,1,10,100],[0,1,10,100])
plt.xlabel('SMA (pix)',size=14)
plt.ylabel('x0 & Y0',size=14)
plt.legend()


"""
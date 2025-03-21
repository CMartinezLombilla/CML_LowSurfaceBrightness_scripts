import numpy as np
import os
import pickle

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.close('all')
plt.ioff()

plt.rcParams["font.family"] = "serif"
plt.rc('text', usetex=False)

from astropy.io import fits
from astropy.table import Table
from photutils import EllipticalAperture


# =============================================================================
# Import CML
# =============================================================================

import sys
sys.path.append('/Users/felipe/Dropbox/CML/') # path del CML Feli
sys.path.append('/Users/z3530031/Dropbox/CML/') # path del CML Cris

sys.path.append('/Users/felipe/Dropbox/CML/WHT/') # path del fichero de parametros Feli
sys.path.append('/Users/z3530031/Dropbox/CML/WHT/') # path del fichero de parametros Cris

from CML import SecImData
from CML import LogNorm
from CML import Counts2Sb
from CML import save_fits_image

# Division entera en "base"
def baseround(x, base=10):
    return base * round(x/base)


# =============================================================================
# Parametros iniciales del grupo/imagen/IGL/aperturas
# =============================================================================
groupID = '400138'
members_IGL = ['BCG', '1660730', '1660646', '1660615']
members_core = ['BCG', '1660730', '1660646', '1660615']
members = ['BCG','1660545', '1660615', '1660646', '1660730', '2301069', 'IGL']

bands = ['g', 'r', 'i']
maxsma = {'1660545':45, 'Core':500, 'IGL':500, '2301069':40}
extra_pa = np.deg2rad(45)
extra_eps = 0.3

# =============================================================================
# Prepare diccionaries
# =============================================================================
table_headers = ['Band', 'TotalFlux', 'IGL_sum_mod', 'IGL_sum_SB26', 'IGL_sum_SB26.5', 'IGL_sum_SB27', 'IGLfrac_2DMod', 'IGLfrac_SB26', 'IGLfrac_SB26.5' , 'IGLfrac_SB27']

TotalFlux = {'Name':'TotalFlux', 'g':0, 'r':0, 'i':0}
IGL_sum_2Dmod = {'Name':'IGL_sum_2Dmod','g':0, 'r':0, 'i':0}
IGLfrac_2DMod = {'Name':'IGLfrac_2DMod','g':0, 'r':0, 'i':0}

IGL_sum_SB = {'Name':'IGL_sum_SB','g':[], 'r':[], 'i':[]}
IGLfrac_SB  = {'Name':'IGLfrac_SB','g':0, 'r':0, 'i':0}
counts_imlim = {'g':0, 'r':0, 'i':0}
counts_IGLlim = {'g':[], 'r':[], 'i':[]}
sb_IGLlim = {'g':[], 'r':[], 'i':[]}

# =============================================================================
# Paths
# =============================================================================
results_path = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis-HSC/results/'+groupID+'/'


# =============================================================================
# Creamos un directorio para guardar los distintos productos del programa
# =============================================================================

try:
    # Create target Directory
    os.mkdir(results_path+'IGL/IGLfraction')
except:
    pass


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

# Parametros del instrumento y la imagen
zp = params.zp
pixscale = params.pixelscale
SkyRMS = params.SkyRMS
HSC_sb_lim = params.sb_limit

# A y Correccion K
k_corrs = params.k_corrs
extinctions = params.extinctions
SkyRMS = params.SkyRMS

# sb dimming: -2.5*np.log10(intens) - 10*np.log10(1+z) -- pq el dimming va como (1+z)**-4 -- Calvi+2014  doi: 10.1088/0004-637X/796/2/102
dimming = 10*np.log10(1+z) 

#Posiciones centrales de cada miembro del grupo, de ls IGL y del core
Xcoords = []
Ycoords = []
for member in members:
    Xcoords.append(int(round(params.Xc[member],0)))
    Ycoords.append(int(round(params.Yc[member],0)))

Xigl = []
Yigl = []
for memberigl in members_IGL:
    Xigl.append(int(round(params.Xc[memberigl],0)))
    Yigl.append(int(round(params.Yc[memberigl],0)))

 
Xcore = []
Ycore = []
for memberc in members_core:
    Xcore.append(int(round(params.Xc[memberc],0)))
    Ycore.append(int(round(params.Yc[memberc],0)))

    
os.system('rm params.py')



# =============================================================================
# Load the 3 masks of the HSC detector
# =============================================================================

HSC_g = fits.getdata(results_path+'masks/intermediate/detectorHSC-g.fits')
HSC_r = fits.getdata(results_path+'masks/intermediate/detectorHSC-r.fits')
HSC_i = fits.getdata(results_path+'masks/intermediate/detectorHSC-i.fits')

HSC_mask = HSC_g + HSC_r + HSC_i



# =============================================================================
# 0 - Total masss of the system
# =============================================================================

# =============================================================================
# Find the center of masss of the core region
# =============================================================================
Gmasked = fits.getdata(results_path+'masks/GroupMasked-i.fits')
weights_list = []

for ii in np.arange(len(members_core)):
    weights_list.append(Gmasked[Ycore[ii], Xcore[ii]])

weights = np.array(weights_list)

com_x = sum(Xcore * weights) / sum(weights)
com_y = sum(Ycore * weights) / sum(weights)

# Center of mass of the system (use it for the extraction of the IGL)
Xc_core = com_x 
Yc_core = com_y


for ind, band in enumerate(bands):
    
    # SB_lim of the images to counts
    counts_imlim[band] = 10**(-0.4*(HSC_sb_lim[band] - zp - 5.*np.log10(pixscale)))

 
    # =============================================================================
    # Read unmasked data, group mask, and add HSC detector mask
    # =============================================================================
    data_nomask = fits.getdata(results_path+'PSF/PSFcorr-'+band+'-SkySubDist.fits')

    data = fits.getdata(results_path+'masks/GroupMasked-'+band+'.fits')
    
    # Nos aseguramos que estan enmascarados con las 3 mascaras del detector
    data[HSC_mask>0] = np.nan  
    
    mask = np.zeros_like(data)
    mask[np.isnan(data)] = np.nan
    
    # =============================================================================
    # Para la elipse del core, usamos los parametros del ajuste original en BCGsky [+ un extra para ajustar mejor segun lo que necesite)]
    # =============================================================================
    ell_sky = pickle.load(open(results_path+'BCGsky/main_model.p', 'rb' ) )
    
    eps = np.median(ell_sky.isolist.eps)+extra_eps
    pa = np.median(ell_sky.isolist.pa)+extra_pa
    a = maxsma['Core'] 
    b = (1-eps)*a
    
    # Get appertures and sum all pixels inside the ellipse    
    # Core
    Core_ap = EllipticalAperture((Xc_core, Yc_core), a=a, b=b,theta=pa)
    Core_sec = SecImData(Core_ap, data, method='center', subpixels=None)
    
    Core_sum = np.nansum(Core_sec[Core_sec>counts_imlim[band]])
    
    # =============================================================================
    # Member 1660545
    # =============================================================================
    Member = '1660545'
    mem_in = members.index(Member)
    
    eps = 0
    pa = 0
    a = maxsma[Member]
    b = (1-eps)*a
    
    Mem545_ap = EllipticalAperture((Xcoords[mem_in], Ycoords[mem_in]), a=a, b=b,theta=pa)
    Mem545_sec = SecImData(Mem545_ap, data, method='center', subpixels=None)
    Mem545_sum = np.nansum(Mem545_sec[Mem545_sec>counts_imlim[band]])

    
    # =============================================================================
    # Member 2301069
    # =============================================================================
    Member = '2301069'
    mem_in = members.index(Member)
    
    eps = 0
    pa = 0
    a = maxsma[Member]
    b = (1-eps)*a
    
    Mem069_ap = EllipticalAperture((Xcoords[mem_in], Ycoords[mem_in]), a=a, b=b,theta=pa)
    Mem069_sec = SecImData(Mem069_ap, data, method='center', subpixels=None)
    Mem069_sum = np.nansum(Mem069_sec[Mem069_sec>counts_imlim[band]])
    
    # =============================================================================
    # Total Flux of the group in counts
    # =============================================================================
        
    TotalFlux[band] = Core_sum + Mem545_sum + Mem069_sum
    #print('Total Flux of the group in counts ('+band+'-band): '+str(TotalFlux))
    
    
    
    
    # =============================================================================
    # Figuras de comprobacion        
    # =============================================================================
    # =============================================================================
    # PDF para las imagenes de control
    # =============================================================================
    pdf = PdfPages(results_path+'IGL/IGLfraction/All_appertures.pdf')
    
    #==============================================================================
    # Comprobacion: plot de las aperturas
    #==============================================================================
    plt.figure(ind,figsize=(8,8))
    plt.title(band+'-band; Group ID: '+groupID,size=18)
    
    # Datos sin enmascarar
    plt.imshow(data_nomask, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=10))
    Core_ap.plot(color='white', ls=':')
    Mem545_ap.plot(color='white', ls=':')
    Mem069_ap.plot(color='white', ls=':')
    
    # =============================================================================
    # Sombreamos las zonas enmascaradas
    # =============================================================================
    from matplotlib import colors
    cmap = colors.ListedColormap(['k'])
    
    shadow = np.nan_to_num(mask,nan=1)
    shadow[shadow==0]=np.nan
    
    plt.imshow(shadow,cmap=cmap, interpolation='nearest', origin='lower',alpha=0.5)
    plt.xlim(np.shape(data)[0]/2-np.shape(data)[0]/3, np.shape(data)[0]/2+np.shape(data)[0]/3)
    plt.ylim(np.shape(data)[0]/2-np.shape(data)[0]/3, np.shape(data)[0]/2+np.shape(data)[0]/3)
 
    

    plt.show()

    for fig in range(0, plt.figure().number): ## will open an empty extra figure :(
        pdf.savefig( fig )
pdf.close()
plt.close('all')
   

    






# =============================================================================
# 1 - IGL fraction from 2D models
# =============================================================================

# =============================================================================
# Find the center of masss of the IGL region
# =============================================================================
Gmasked = fits.getdata(results_path+'masks/GroupMasked-i.fits')
weights_list = []

for ii in np.arange(len(members_IGL)):
    weights_list.append(Gmasked[Yigl[ii], Xigl[ii]])

weights = np.array(weights_list)

com_x = sum(Xigl * weights) / sum(weights)
com_y = sum(Yigl * weights) / sum(weights)

# Center of mass of the system (use it for the extraction of the IGL)
Xc_IGL = com_x 
Yc_IGL = com_y


for ind, band in enumerate(bands):
    
    # =============================================================================
    # Leemos los datos enmascarados y creamos una mascara
    # =============================================================================
    data_nomask = fits.getdata(results_path+'PSF/PSFcorr-'+band+'-SkySubDist.fits')
    data_IGL = fits.getdata(results_path+'masks/Dist_IGLMasked-'+band+'.fits')
   
    # Nos aseguramos que estan enmascarados con las 3 mascaras del detector
    data_IGL[HSC_mask>0] = np.nan  
    
    mask = np.zeros_like(data_IGL)
    mask[np.isnan(data_IGL)] = np.nan
    
    # =============================================================================
    # Para la elipse de la IGL, usamos los parametros del ajuste original en BCGsky
    # =============================================================================
    ell_sky = pickle.load(open(results_path+'BCGsky/main_model.p', 'rb' ) )
    
    eps = np.median(ell_sky.isolist.eps)+extra_eps
    pa = np.median(ell_sky.isolist.pa)+extra_pa
    a = maxsma['IGL']
    b = (1-eps)*a
    
    # Obtenemos las aperturas y sumamos los pixels de adentro      
    # Core
    IGL_ap = EllipticalAperture((Xc_IGL, Yc_IGL), a=a, b=b,theta=pa)
    IGL_sec = SecImData(IGL_ap, data_IGL, method='center', subpixels=None)
    
    IGL_sum_2Dmod[band] = np.nansum(IGL_sec[IGL_sec>counts_imlim[band]])
    
    # =============================================================================
    # Total Flux of the IGL in counts - 2D models
    # =============================================================================
    #print('Total Flux of the IGL in counts ('+band+'-band): '+str(IGL_sum_mod))
    
    
    # =============================================================================
    # Fraction of IGL - 2D Models
    # =============================================================================
    
    IGLfrac_2DMod[band] = IGL_sum_2Dmod[band]/TotalFlux[band]
    print('***************************************************************')
    print('IGL fraction from 2D Models ('+band+'-band): '+str(IGLfrac_2DMod[band]))
    print('***************************************************************')

    
    # =============================================================================
    # Figuras de comprobacion        
    # =============================================================================
    # =============================================================================
    # PDF para las imagenes de control
    # =============================================================================
    pdf = PdfPages(results_path+'IGL/IGLfraction/IGL-2DMod_ap'+str(maxsma['IGL'])+'pix.pdf')
    
    #==============================================================================
    # Comprobacion: plot de las aperturas
    #==============================================================================
    plt.figure(ind,figsize=(8,8))
    plt.title('Group ID: '+groupID+'; '+band+'-band',size=18)
    
    # Datos sin enmascarar
    plt.imshow(data_nomask, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=10))
    IGL_ap.plot(color='white', ls=':')
    
    # =============================================================================
    # Sombreamos las zonas enmascaradas
    # =============================================================================
    cmap = colors.ListedColormap(['k'])
    
    shadow = np.nan_to_num(mask,nan=1)
    shadow[shadow==0]=np.nan
    
    plt.imshow(shadow,cmap=cmap, interpolation='nearest', origin='lower',alpha=0.5)
    plt.xlim(np.shape(data_IGL)[0]/2-np.shape(data_IGL)[0]/3, np.shape(data_IGL)[0]/2+np.shape(data_IGL)[0]/3)
    plt.ylim(np.shape(data_IGL)[0]/2-np.shape(data_IGL)[0]/3, np.shape(data_IGL)[0]/2+np.shape(data_IGL)[0]/3)
 
    

    plt.show()

    for fig in range(0, plt.figure().number): ## will open an empty extra figure :(
        pdf.savefig( fig )
pdf.close()
plt.close('all')


###############################################################################
###############################################################################



#bands = ['g']


# =============================================================================
# 2 - IGL fraction from sb limits
# =============================================================================

# SB LIMITS (down_sb_lim = HSC_sb_lim)
sb_IGLlim_rf = [24, 24.5, 25, 25.5, 26, 26.5, 27, 27.5, 28]
#bands = ['i']

for ind, band in enumerate(bands):
 
    # =============================================================================
    # Leemos los datos enmascarados y creamos una mascara
    # =============================================================================
    data_nomask = fits.getdata(results_path+'PSF/PSFcorr-'+band+'-SkySubDist.fits')
    data_IGL = fits.getdata(results_path+'masks/CoreMasked-'+band+'.fits')
    # Nos aseguramos que estan enmascarados con las 3 mascaras del detector
    data_IGL[HSC_mask>0] = np.nan  
    
    mask = np.zeros_like(data_IGL)
    mask[np.isnan(data_IGL)] = np.nan
    
    
    
    # =============================================================================
    # =============================================================================
    # METHOD 2: surface brightness limits (convert to counts of the sb_lims values)
    # =============================================================================
    # =============================================================================
    IGL_sum_SBi = []
    IGLfrac_SBi = []
       
    # masking below image lim (HSC_sb_lim)
    data_down = np.copy(data_IGL)
    data_down[data_IGL<counts_imlim[band]] = np.nan
    
        
        
    # Masking above upper lim (IGL lim) for each sb_uplim
    for ii, sblim in enumerate(sb_IGLlim_rf):
        
        sb_IGLlim[band].append(sblim)# + k_corrs[band] + dimming + extinctions[band])
        
        counts_IGLlim[band].append(10**(-0.4*(sb_IGLlim[band][ii] - zp - 5.*np.log10(pixscale))))
        
        data_IGL_ii = np.copy(data_down)
        
        data_IGL_ii[data_IGL>counts_IGLlim[band][ii]] = np.nan
        
        
        # Get appertures and sum all pixels inside the apperture of the core region    
        SBlim_sec = SecImData(Core_ap, data_IGL_ii, method='center', subpixels=None)
        IGL_sumii = np.nansum(SBlim_sec)
        IGL_sum_SBi.append(IGL_sumii)
    
        # =============================================================================
        # Fraction of IGL - at each sb limit
        # =============================================================================     
        IGLfrac_SBi.append(IGL_sumii/TotalFlux[band])

        # Save fits image with the fraction of light measured
        save_fits_image(data_IGL_ii, fits.getheader(results_path+'PSF/PSFcorr-'+band+'-SkySubDist.fits'), results_path+'IGL/IGLfraction/' , 'IGLfrac-SBlim-'+str(sblim)+'-'+band+'.fits')


        '''
        # =============================================================================
        # Figuras de comprobacion        
        # =============================================================================
        
        pdf = PdfPages(results_path+'IGL/IGLfraction/IGL-SBlim-'+str(sblim)+'-'+band+'-Images.pdf')

        #==============================================================================
        # Plot de las aperturas
        #==============================================================================
        plt.figure(1, figsize=(8,8))
        plt.title('Group ID: '+groupID+'; '+band+'-band; SBlim '+str(sblim),size=18)
        
        # Datos sin enmascarar
        plt.imshow(data_IGL_ii, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=3))        
        # =============================================================================
        # Sombreamos las zonas enmascaradas
        # =============================================================================
        cmap = colors.ListedColormap(['k'])
        
        #shadow = np.zeros_like(data_counts)
        #shadow[np.isnan(data_counts)] = 1
        #shadow[shadow==0]=np.nan
        
        #plt.imshow(shadow,cmap=cmap, interpolation='nearest', origin='lower',alpha=0.5)
        plt.xlim(np.shape(data_IGL)[0]/2-np.shape(data_IGL)[0]/3, np.shape(data_IGL)[0]/2+np.shape(data_IGL)[0]/3)
        plt.ylim(np.shape(data_IGL)[0]/2-np.shape(data_IGL)[0]/3, np.shape(data_IGL)[0]/2+np.shape(data_IGL)[0]/3)
        
        pdf.savefig(1)
        pdf.close()
        plt.close('all')
        '''

    # guardamos los valores finales en cada banda
    IGL_sum_SB[band] = IGL_sum_SBi
    IGLfrac_SB[band] = IGLfrac_SBi  
     

  

'''
    # =============================================================================
    # =============================================================================
    # METHOD 1: counts to surface brightness over the whole images
    # =============================================================================
    # =============================================================================
    IGL_sum_SBi = []
    IGLfrac_SBi = []

    c2sb = Counts2Sb(zp, pixscale, A=extinctions[band], dimming=dimming, k_corr=k_corrs[band])
    #c2sb = Counts2Sb(zp, pixscale, A=0, dimming=0, k_corr=0)

    data_IGL_sb = c2sb.do(data_IGL)
    
    # masking below down lim (HSC_sb_lim)
    data_IGL_sb[data_IGL_sb>HSC_sb_lim[band]] = np.nan
    
    # =============================================================================
    # PDF para las imagenes de control
    # =============================================================================
        
    pdf = PdfPages(results_path+'IGL/IGLfraction/IGL-SBlim-'+band+'-Images.pdf')
        
        
    # Masking above upper lim for each sb_uplim
    for ii, sblim in enumerate(sb_uplim):
        
        data_IGL_sbii = np.copy(data_IGL_sb)
        data_IGL_sbii[data_IGL_sb<sblim] = np.nan
        
        # Sumamos los pixels que nos quedan (en cuentas!)
        data_counts = np.copy(data_IGL)
        data_counts[np.isnan(data_IGL_sbii)] = np.nan
        
        IGL_sumii = np.nansum(data_counts)
        IGL_sum_SBi.append(IGL_sumii)
    
        # =============================================================================
        # Fraction of IGL - at each sb limit
        # =============================================================================     
        IGLfrac_SBi.append(IGL_sumii/TotalFlux[band])



        # =============================================================================
        # Figuras de comprobacion        
        # =============================================================================
        #==============================================================================
        # Plot de las aperturas
        #==============================================================================
        plt.figure(1, figsize=(8,8))
        plt.title('Group ID: '+groupID+'; '+band+'-band; SBlim '+str(sblim),size=18)
        
        # Datos sin enmascarar
        plt.imshow(data_counts, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=3))        
        # =============================================================================
        # Sombreamos las zonas enmascaradas
        # =============================================================================
        cmap = colors.ListedColormap(['k'])
        
        #shadow = np.zeros_like(data_counts)
        #shadow[np.isnan(data_counts)] = 1
        #shadow[shadow==0]=np.nan
        
        #plt.imshow(shadow,cmap=cmap, interpolation='nearest', origin='lower',alpha=0.5)
        plt.xlim(np.shape(data_IGL)[0]/2-np.shape(data_IGL)[0]/3, np.shape(data_IGL)[0]/2+np.shape(data_IGL)[0]/3)
        plt.ylim(np.shape(data_IGL)[0]/2-np.shape(data_IGL)[0]/3, np.shape(data_IGL)[0]/2+np.shape(data_IGL)[0]/3)
        
        pdf.savefig(1)
        
        plt.close('all')
     
    pdf.close()
    
    # guardamos los valores finales en cada banda
    IGL_sum_SB[band] = IGL_sum_SBi
    IGLfrac_SB[band] = IGLfrac_SBi   
     
'''    




print('***************************************************************')
print('IGL fraction from SB limits: '+str(IGLfrac_SB))
print('***************************************************************')



          
# =============================================================================
#  Save a table with all the info  
# =============================================================================

table_vec = [TotalFlux, IGL_sum_2Dmod, IGLfrac_2DMod]

for ii, sblim in enumerate(sb_IGLlim_rf):
    frac_sb = {'Name':'IGLfrac_SB '+str(sblim),'g':IGLfrac_SB['g'][ii], 'r':IGLfrac_SB['r'][ii], 'i':IGLfrac_SB['i'][ii]}
    table_vec.append(frac_sb)
    
        
table = Table(rows=table_vec, meta={'IGL amount and fraction': 'table info'})
table.write(results_path+'IGL/IGLfraction/IGLfrac_table.fits', overwrite=True)


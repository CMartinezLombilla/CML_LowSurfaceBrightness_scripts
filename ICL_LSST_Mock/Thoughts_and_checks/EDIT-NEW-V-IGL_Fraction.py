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


from CML import SecImData
from CML import LogNorm
from CML import save_fits_image
from CML import Center


# =============================================================================
# Initial params
# =============================================================================
sim = 'Horizon_AGN'

bands = ['r']
BCG_coords = [10266, 10284] # written as [Y_ds9, X_ds9]

pixscale = 0.2
lsst_sb_lims = {'u':29.4, 'g':30.3, 'r':30.3, 'i':29.7, 'z':28.9, 'y':28.1}


# =============================================================================
# Define diccionaries
# =============================================================================
 
cluster_total_flux = {'Name':'TotalFlux', 'g':0, 'r':0, 'i':0}
ICL_sum_2Dmod = {'Name':'ICL_sum_2Dmod','g':0, 'r':0, 'i':0}
ICLfrac_2DMod = {'Name':'ICLfrac_2DMod','g':0, 'r':0, 'i':0}

icl_sum_sb = {'Name':'ICL_sum_SB','g':[], 'r':[], 'i':[]}
icl_frac_sb  = {'Name':'ICLfrac_SB','g':0, 'r':0, 'i':0}
counts_imlim = {'g':0, 'r':0, 'i':0}
counts_icl_lim = {'g':[], 'r':[], 'i':[]}


# =============================================================================
# Paths 
# =============================================================================
data_dir = '/Users/z3530031/unsw/Work/Projects/LSST/Chall3-ICL/data/'+sim+'/'
res_dir = '/Users/z3530031/unsw/Work/Projects/LSST/Chall3-ICL/results/'+sim+'/'

image_file = '00761_0000048_0.05_xy_r_4Mpc.fits'
#image_files = glob.glob(image_dir+'*.fits') # Load a list 


# =============================================================================
# Create directory to save results
# =============================================================================
try:
    # Create target Directory
    os.mkdir(res_dir+'ICL_fraction')
except:
    pass




# =============================================================================
# 0 - Total masss of the system
# =============================================================================


for ind, band in enumerate(bands):
    
    # LSST SB_lim of the images to counts
    counts_imlim[band] = 10**(-0.4*(lsst_sb_lims[band] - 5.*np.log10(pixscale)))
 
    # Load images
    im, hdr = fits.getdata(data_dir+image_file, header = True)
        
    
    # =============================================================================
    # Find the center of the BCG in a small region of 10 pixels
    # =============================================================================
    
    center=Center(im, BCG_coords[0], BCG_coords[1], 10)
    
    # Fixed centres
    Xc=center.x
    Yc=center.y

    # =============================================================================
    # Total flux in the image in counts (total cluster flux)
    # =============================================================================
    cluster_total_flux[band] = np.nansum(im[im>counts_imlim[band]])
    print('Total Flux of the group in counts ('+band+'-band): '+str(cluster_total_flux[band]))
    
    
    
'''   
    # =============================================================================
    # Figuras de comprobacion        
    # =============================================================================

    pdf = PdfPages(results_path+'IGL/IGLfraction/All_appertures.pdf')
    
    #==============================================================================
    # Comprobacion: plot de las aperturas
    #==============================================================================
    plt.figure(ind,figsize=(8,8))
    plt.title(band+'-band; sim: '+sim,size=18)
    
    # Datos sin enmascarar
    plt.imshow(data_nomask, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=10))

    # Sombreamos las zonas enmascaradas
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
   
'''
    




'''

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



 
'''


# =============================================================================
# 2 - ICL fraction from sb limits
# =============================================================================

# SB LIMITS (down_sb_lim = HSC_sb_lim)
sb_icl_lims = [26, 26.5, 27]
#bands = ['i']

for ind, band in enumerate(bands):
 
    # Load images
    im, hdr = fits.getdata(data_dir+image_file, header = True)
    
    #data_nomask = fits.getdata(results_path+'PSF/PSFcorr-'+band+'-SkySubDist.fits')
    #data_IGL = fits.getdata(results_path+'masks/CoreMasked-'+band+'.fits')
    
    #mask = np.zeros_like(data_IGL)
    #mask[np.isnan(data_IGL)] = np.nan
    
    
    
    # =============================================================================
    # METHOD: convert to counts the the sb_lims values
    # =============================================================================
    icl_sum_sbi = []
    icl_frac_sbi = []
       
    # masking below image lim (HSC_sb_lim)
    data_down = np.copy(im)
    data_down[im<counts_imlim[band]] = np.nan 
        
    
    for ii, sblim in enumerate(sb_icl_lims):
        
        # Masking above upper lim (ICL lim) for each sb_uplim
        counts_icl_lim[band].append(10**(-0.4*(sblim - 5.*np.log10(pixscale))))
        data_icl_ii = np.copy(data_down)
        data_icl_ii[im>counts_icl_lim[band][ii]] = np.nan
                
        # Sum all the unmasked pixels    
        icl_sumii = np.nansum(data_icl_ii)
        icl_sum_sbi.append(icl_sumii)
    
        # =============================================================================
        # Fraction of ICL at each sb limit
        # =============================================================================     
        icl_frac_sbi.append(icl_sumii/cluster_total_flux[band])

        # Save fits image with the fraction of light measured
        save_fits_image(data_icl_ii, hdr, res_dir+'ICL_fraction/', 'SBlim-'+str(sblim)+'-'+image_file)

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
    icl_sum_sb[band] = icl_sum_sbi
    icl_frac_sb[band] = icl_frac_sbi  
     



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

'''
# =============================================================================
# Imports...
# =============================================================================
import numpy as np
import os
from astropy.table import Table
from astropy.wcs import WCS
from astropy.io import fits
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from adjustText import adjust_text # conda install -c conda-forge adjusttext 

# =============================================================================
# Import CML
# =============================================================================
import sys
sys.path.append('/Users/felipe/Dropbox/CML/') # path del CML Feli
sys.path.append('/Users/z3530031/Dropbox/CML/') # path del CML Cris

from CML import LogNorm



# =============================================================================
# Initial params and paths
# =============================================================================


data_dir = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis/data/'
results_dir = '/Users/z3530031/unsw/Work/Projects/ICL_HSC-LSST/Analysis/results/'


# Tables with data I'll use
tables_dir = data_dir+'tables/'

gal_table_path = tables_dir+'G3CGalv10.fits'  
group_table_path = tables_dir+'G3CFoFGroupv10.fits'
extinc_table_path = tables_dir+'UDeep_G02GalacticExtinctionv01.fits'  
kcorr_table_path = tables_dir+'UDeep_G02kcorr_z00_v01.fits'  




# =============================================================================
# Create a directory for the group and some subfolders
# =============================================================================


#GroupID = 400016
groupID_list_yes = [400016,
                    400020,
                    400026,
                    400030,
                    400032,
                    400034,
                    400039,
                    400043,
                    400049,
                    400052,
                    400054,
                    400069,
                    400070,
                    400073,
                    400081,
                    400102,
                    400103,
                    400109,
                    400115,
                    400119,
                    400125,
                    400128,
                    400138,
                    400145,
                    400146,
                    400153,
                    400162,
                    400164,
                    400165,
                    400170,
                    400172,
                    400173,
                    400217,
                    400218,
                    400219,
                    400221,
                    400230,
                    400231,
                    400241,
                    400255,
                    400258,
                    400275,
                    400282,
                    400298,
                    400365,
                    400367,
                    400386,
                    400393]  


for GroupID in groupID_list_yes:

    try:
        # Create target Directory
        os.mkdir(results_dir+str(GroupID))
    except FileExistsError:
        print ('# ==============================')
        print ('# Group directory already exist')
        print ('# ==============================')
    
    try:
        # Create target Directory
        os.mkdir(results_dir+str(GroupID)+'/masks')
    except FileExistsError:
        print ('# ==============================')
        print ('# masks directory already exist')
        print ('# ==============================')
    
    try:
        # Create target Directory
        os.mkdir(results_dir+str(GroupID)+'/masks/intermediate')
    except FileExistsError:
        print ('# ==============================')
        print ('# masks/intermediate directory already exist')
        print ('# ==============================')
    
               
    # =============================================================================
    # Paths for the original images
    # =============================================================================
    
    images_dir = data_dir+'images/'
    data_path = images_dir+str(GroupID)+'-gri.fits'
    
    # =============================================================================
    # Instrument params
    # =============================================================================
    pixelscale = 0.168 # arcsec/pixel
    zp = 27 #mag for all HSC PDR2 griyz bands  -- LSST STACK: zp_flux=im.getPhotoCalib().getInstFluxAtZeroMagnitude() ; zp = -2.5*np.log10(zp_flux)   # zero point in mag
    
    # =============================================================================
    # Other params from the tables
    # =============================================================================
    # Parametros del grupo 
    #  GroupID for looking the CATAID en G3CGalv10.fits
    
    # Table filtered with only the member in the group
    gal_table = Table.read(gal_table_path)
    gal_table = gal_table[gal_table['GroupID']==GroupID]
    
    
    # The BCG is the memeber with RankBCG = 1
    bcgID = gal_table[gal_table['RankBCG']==1]['CATAID'][0]
    
    # Use z from the group's table
    group_table = Table.read(group_table_path)
    group_table = group_table[group_table['GroupID']==GroupID]
    
    z = group_table['IterCenZ'][0]
    
    # sb dimming: -2.5*np.log10(intens) - 10*np.log10(1+z) -- pq el dimming va como (1+z)**-4 -- Calvi+2014  doi: 10.1088/0004-637X/796/2/102
    dimming = 10*np.log10(1+z) 
    
    # List group members' CATIDs  
    cata_id = np.sort(gal_table['CATAID'])
    
    
    # =============================================================================
    # Galactic extinctions and k-corrections
    # =============================================================================
    
    # Galactic extinctions for each band
    extinc_table = Table.read(extinc_table_path)
    extinc_table = extinc_table[extinc_table['CATAID']==bcgID]
    
    extinc_g = extinc_table['A_g'][0]
    extinc_r = extinc_table['A_r'][0]
    extinc_i = extinc_table['A_i'][0]
    
    
    # K-corrections for each band
    kcorr_table = Table.read(kcorr_table_path)
    kcorr_table = kcorr_table[kcorr_table['CATAID']==bcgID]
    
    kcorr_g = kcorr_table['KCORR_G'][0]
    kcorr_r = kcorr_table['KCORR_R'][0]
    kcorr_i = kcorr_table['KCORR_I'][0]
    
    
    # =============================================================================
    # Center of each galaxy in the group 
    # =============================================================================
    # Read astrometry
    wcs = WCS(fits.getheader(data_path))
    
    # Create dictionaries with the initial centers
    Xc ={}
    Yc ={}
    
    # The first element is the BCG 
    ra = gal_table[gal_table['CATAID']==bcgID]['RA'].data[0] # in degrees
    dec = gal_table[gal_table['CATAID']==bcgID]['Dec'].data[0] # in degrees
    coord = SkyCoord(ra, dec, unit="deg")
    x, y = coord.to_pixel(wcs)
    # Set the values and keys in the dicctionary for BCG and BCG-sky
    Xc = {'BCG-sky':float(x),'BCG':float(x)} 
    Yc = {'BCG-sky':float(y),'BCG':float(y)}
    
    
    
    # =============================================================================
    # Build the remaining dictionaries
    # =============================================================================
    mask_file = {'BCG-sky':'BCGMask.fits','BCG':'BCGMask.fits'} 
    
    
    for i in cata_id:
        # Asign the key
        if i != bcgID:
            key = str(i)
        
            ra = gal_table[gal_table['CATAID']==i]['RA'].data[0] # in degrees
            dec = gal_table[gal_table['CATAID']==i]['Dec'].data[0] # in degrees
        
            coord = SkyCoord(ra, dec, unit="deg")
            # Proyect the coordinate in the image
            x, y = coord.to_pixel(wcs)
            # Set the value in the positions' diccionary
            Xc[key] = float(x)
            Yc[key] = float(y)
            # Diccionary for the masks
            mask_file[key] = 'GroupMask.fits'
    
    # =============================================================================
    # Diccionaries for the ellipse fitting with some random initial values
    # =============================================================================
    sma0 = mask_file.copy()
    maxsma = mask_file.copy()
    step = mask_file.copy()
    ell = mask_file.copy() 
    theta = mask_file.copy() 
    
    for i in sma0:
        sma0[i] = 5
        maxsma[i] = 200
        step[i] = 0.2
        ell[i] = 0.2 
        theta[i] = 0
    
    
    
    # =============================================================================
    # Plot the positions
    # =============================================================================
    plt.close('all')
    plt.figure(1,figsize=(8,8))
    plt.title('Group ID: '+str(GroupID),size=25)
    
    data = fits.getdata(data_path)
    plt.imshow(data, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=10))
    
    
    
    # Make lists for members' labels and their coords
    label = [*Xc][1:]
    x = []
    y = []
    
    for i in label:    
        # Plot circles in all the positions
        plt.scatter(Xc[i],Yc[i], facecolors='none', s=10,color='w',lw=0.4)
        x.append(Xc[i])
        y.append(Yc[i])
    
    texts = [plt.text(x[i], y[i],label[i], va='center',ha='center',color='w') for i in range(len(label))]
    adjust_text(texts)
    
    plt.xticks([],[])
    plt.yticks([],[])
    plt.tight_layout()
    
    plt.savefig(results_dir+str(GroupID)+'/GroupMap.pdf',dpi=300)
    
    
    
    
    # =============================================================================
    # TO DO: Create a template
    # =============================================================================
    
    # =============================================================================
    # Read the template
    # =============================================================================
    with open(data_dir+'params_template.py', 'r') as file :
      filedata = file.read()
    # Params of the grup
    filedata = filedata.replace('ID_GRUPO', str(GroupID))
    filedata = filedata.replace('ID_BCG', str(bcgID))
    filedata = filedata.replace('REDSHIFT', str(z))
    # Params of the instrument
    filedata = filedata.replace('ESCALA_DEL_PIXEL', str(pixelscale))
    filedata = filedata.replace('PO_POINT', str(pixelscale))
    # From the tables...
    filedata = filedata.replace('EXTINCTIONS', str({'g':extinc_g, 'r':extinc_r, 'i':extinc_i}))
    filedata = filedata.replace('K-CORRECTIONS', str({'g':kcorr_g, 'r':kcorr_r, 'i':kcorr_i}))
    # Params for each galaxy for the ellipse fitting
    filedata = filedata.replace('DICT_XC', str(Xc).replace(',',',\n'))
    filedata = filedata.replace('DICT_YC', str(Yc).replace(',',',\n'))
    filedata = filedata.replace('DICT_MASK_FILE', str(mask_file).replace(',',',\n'))
    filedata = filedata.replace('DICT_SMA0', str(sma0).replace(',',',\n'))
    filedata = filedata.replace('DICT_MAXSMA', str(maxsma).replace(',',',\n'))
    filedata = filedata.replace('DICT_STEP', str(step).replace(',',',\n'))
    filedata = filedata.replace('DICT_ELL', str(ell).replace(',',',\n'))
    
    
    filedata = filedata.replace('DICT_THETA', str(theta).replace(',',',\n'))
    
    
    
    # Write parms file
    with open(results_dir+str(GroupID)+'/params.py', 'w') as file:
      file.write(filedata)
    
    












#*****************************************************************************
#                             SB PROFILES FUNCTIONS 
#-----------------------------------------------------------------------------
#*****************************************************************************

# Packages versions
import numpy as np
import matplotlib.pyplot as plt 

# Phoutils
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import EllipticalAperture

# Astropy
from astropy.stats import sigma_clipped_stats        
from astropy.table import Table   
from astropy.io import fits
from astropy.wcs import WCS

# Scipy
from scipy.optimize import minimize_scalar


# Versions
# numpy== 1.19.2
# matplotlib==3.3.2
# astropy==4.1
# photutils==1.0.1
# tkinter == 8.6

# conda install astropy==4.1
# conda install photutils==1.0.1


#==============================================================================
# Ignore Warnings
#==============================================================================
import warnings
warnings.filterwarnings("ignore")

# Avoid warning if the data is masked, this is, an single image with data + NaNs
warnings.warn("Input data contains invalid values (NaNs or infs), which were automatically clipped.") 


# =============================================================================
# Anciliary Funcions
# =============================================================================

class Counts2Sb:
    def __init__(self, zp, pixscale, A=0, dimming=0, k_corr=0):
        '''
        Class to prepare the convertion from counts to surface brightness [mag/arcsec^2]

        Parameters
        ----------
        zp : float
            Zero point [mag]
        A : float, optional
            Galactic absorption or extinction for a given band. The default is 0.
        dimming : float, optional
            Dimming: 10*np.log10(1+z) [mag/arcsec^2]. 
            Because the diming in linear scale is (1+z)**-4 [Calvi+2014, doi: 10.1088/0004-637X/796/2/102]. 
            The default is 0.
        k_corr : float, optional
            K correction [J. Loveday et al. (2012)]. When obs. redshifted galaxies there is an spectral shift, this is a change of the wavelength values as lambda*(1+z) so as wavelenth range. 
            The default is 0.

        Returns
        -------
         Class object.

        '''
        self.zp = zp
        self.A = A
        self.dimming = dimming
        self.k_corr = k_corr
        self.pixscale = pixscale
        
    def do(self, counts):
        '''
        Function to convert 'x' counts to surface brightness [mag/arcsec^2] with the parameters from Counts2Sb.

        Parameters
        ----------
        x : float or numpy array
            Values in counts.

        Returns
        -------
        Float or numpy array in [mag/arcsec^2].

        '''
        self.counts = counts
        #print(counts,- 2.5*np.log10(counts) + self.zp + 5*np.log10(self.pixscale) - self.A - self.dimming - self.k_corr)
        return - 2.5*np.log10(counts) + self.zp + 5*np.log10(self.pixscale) - self.A - self.dimming - self.k_corr 
    
    def do_errors(self, counts, counts_err):
        
        self.counts = counts 
        up = self.counts + counts_err
        down = self.counts - counts_err
        
        # Pasamos a magnitudes
        error_down = -2.5*np.log10(self.counts/up) 
        error_up = -2.5*np.log10(down/self.counts)
        return [error_up, error_down]




def sb2counts(sb, zp, pixscale, A=0, dimming=0, k_corr=0):
    '''
    Function to convert 'x' surface brightness [mag/arcsec^2] to counts
    
    Parameters
    ----------
    sb : float or numpy array
        Values in counts.
    
    Returns
    -------
    Float or numpy array in surface brightness [mag/arcsec^2]
    
    '''

    return 10**(-0.4*(sb - zp - 5.*np.log10(pixscale) + A + dimming + k_corr) )
    
    
        
        
    
def SecImData(ap,image,method,subpixels=None):
    """
    Extracts a section (using any astropy.aperture shape) of a 2D array

    Parameters
    ----------
    ap : astropy.aperture class
        Apperture class to define a geometry.
    image : numpy array
        Data where to extract a section.
    method : interpolation method
         'exact': the exact intersection of the aperture with each pixel is calculated. 
         'center': a pixel is considered to be entirely in or out of the aperture depending on whether its center is in or out of the aperture. 
         'subpixel': pixels are divided into a number of subpixels, which are in or out of the aperture based on their centers. For this method, the number of subpixels needs to be set with the subpixels keyword.
         'center' and 'subpixel', are faster, but with the expense of less precision.
    subpixels : int, optional
        If 'subpixel' method is used. The default is None.

    Returns
    -------
    sec_data : 1D numpy array
        1D data within the aperture weighted according to the method, with no NaNs, ready for any statistic analysis.

    """
    mask=ap.to_mask(method=method,subpixels=subpixels)
    mask_data=mask.data
    sec=mask.cutout(image)
    sec_weight=sec*mask_data
    sec_data = sec_weight[mask_data>0]
    # Quitamos los valores NaN
    sec_data = sec_data[np.isfinite(sec_data)]
    
    return sec_data


def com(data, Xcoords, Ycoords):
    """
    Function to obtain the center of 'light' (mass) for a set of positions in 
    an astronomical image (or any array). 
    - i.e., to get the center of 'light' (mass) of a system of 3 galaxies by 
    giving their central coords

    Parameters
    ----------
    data : np.array
        Image with the sources
    Xcoords : list
        X coordinates of the sources as for ds9.
    Ycoords : list 
        Y coordinates of the sources as for ds9.

    Returns
    -------
    list
        [x, y] coordinates of the center of 'light' (mass) as for Python.

    """
    
    weights_list = []

    for ii in np.arange(len(Xcoords)):
        weights_list.append(data[Ycoords[ii], Xcoords[ii]])

    weights = np.array(weights_list)
    
    # Center of light of the system 
    com_x = sum(Xcoords * weights) / sum(weights)
    com_y = sum(Ycoords * weights) / sum(weights)

    return [com_x, com_y]




class Center:
    # Ajuste para buscar el centro
    def __init__(self,image,Xc,Yc,d):
        
        self.image=image
        self.Xc=Xc
        self.Yc=Yc
        self.d=d

        xmas=np.int(Xc+d)
        xmenos=np.int(Xc-d)  
        ymas=np.int(Yc+d)
        ymenos=np.int(Yc-d)
        
        #Construimos una seccion para calcular el centroide
        self.sec=image[ymenos:ymas,xmenos:xmas]
        
        from photutils import centroid_2dg # centroid_1dg, centroid_2dg centroid_com
        
        #Las posiciones de los centroides ajustados
        x_s, y_s =centroid_2dg(self.sec)
        
        self.x=x_s+self.Xc-d
        self.y=y_s+self.Yc-d

                
    def find(self):              
        return(self.x,self.y)
    

        
    def plot(self):
        p_x=np.sum(self.sec,0)
        p_y=np.sum(self.sec,1)
        ly=np.arange(len(p_y))+self.Yc-self.d
        lx=np.arange(len(p_x))+self.Xc-self.d
        
        plt.close('all')
        plt.figure(1,figsize=(8,4))

        # Perfil en eje X
        plt.subplot(121)
        plt.title('X direction')
        plt.step(lx,p_x)
        plt.axvline(self.x,color='r',ls='dashed')
        plt.xlabel('pix')
        plt.ylabel('Counts')
        
        # Perfil en eje Y
        plt.subplot(122)
        plt.title('Y direction')
        plt.step(ly,p_y)
        plt.axvline(self.y,color='r',ls='dashed')
        plt.xlabel('pix')
        plt.yticks([],[])
        
        plt.show()
        


#********************************************************************************
# RADIAL CENTRAL PROFILE 
#********************************************************************************


#==============================================================================
# Perfiles radiales. Funcion raiz
#==============================================================================


class CentralProfile:
    def __init__(self,image,Xc,Yc,nbins,npix,height,zp,pixscale,A=0,rms_sky=0,orientation='horizontal',delta=0):
        """
        Main class definition to extract a surface brightness profile using rectangular apertures over an image. Here the apertures are built.     
        
        Parameters
        ----------
        image : 2D numpy array
            Data where to extract the profile.
        Xc : float
            X center coordinate in pixels.
        Yc : float
            Y center coordinate in pixels.
        nbins : int
            Number of bins of the profile.
        npix : int
            The full width of the profile in pixels. Number of pixels up to where extract the profile.
        height : float
            The full height of the rectangles in pixels.
        zp : float
            Zero point [mag/arcsec^2]. The default is 0.
        A : float
            Galactic absorption or extinction for a given band. The default is 0.
        rms_sky : float
            RMS of the sky to include as a source of error. The default is 0.
        orientation : str, optional
            Orientation of the profile. 
            The default is 'horizontal'. 'horizontal' is along the galaxy major axis; 'vertical' if along the minor axis.
        delta : float, optional
            Offset in pixels from the galaxy center where to extract the profile. The default is 0, i.e., along the main axis, either 'horizontal' or 'vertical'.

        Returns
        -------
        Class object with the apertures already defined.

        """
        # Geometric parameters
        self.image=image
        self.Xc = Xc
        self.Yc = Yc
        self.nbins = nbins
        self.npix = npix
        self.height = height
        self.delta=delta
        self.orientation=orientation
        
        # Data parametes
        self.zp = zp
        self.A = A
        self.rms_sky = rms_sky
        self.pixscale = pixscale
        
        #=======================================================================================================
        #  Apertures definition
        #=======================================================================================================
        import numpy as np
        l=np.log10(self.npix)
        step=np.logspace(0,l,self.nbins)

        # Diccionario para todos los parametros segun la apertura sea vertical u horizontal.
        self.param={'horizontal':{'w' : np.diff(step), 'h' : np.ones(len(np.diff(step)))*self.height, 'x' : step[1:]-(np.diff(step))/2-step[0]-np.diff(step)[0]/2, 'y' : np.zeros(len(np.diff(step)))+self.delta, 'signo' : (-1,1)},            
               'vertical':{'w' : np.ones(len(np.diff(step)))*self.height,'h' : np.diff(step), 'x' : np.zeros(len(np.diff(step)))+self.delta, 'y' : step[1:]-(np.diff(step))/2-step[0]-np.diff(step)[0]/2, 'signo' : (1,-1)}}

        w = self.param[orientation]['w']
        h = self.param[orientation]['h']
        x = self.param[orientation]['x']
        y = self.param[orientation]['y']
        self.signo = self.param[orientation]['signo']

        #==============================================================================
        # Apertures built
        #==============================================================================
        from photutils import RectangularAperture

        self.ap_c=RectangularAperture([x[0]+self.Xc,y[0]+self.Yc], w=w[0], h=h[0],theta=0) # Apertura central
        self.ap_r=[] # Aperturas hacia la derecha
        self.ap_l=[] # Aperuras hacia la izquierda
        
                
        for x0, y0, w0, h0 in zip(x[1:], y[1:], w[1:], h[1:]):
            self.ap_r.append(RectangularAperture((x0+self.Xc,y0+self.Yc), w=w0, h=h0,theta=0))
            self.ap_l.append(RectangularAperture((self.signo[0]*x0+self.Xc,self.signo[1]*y0+self.Yc), w=w0, h=h0,theta=0))

        self.w = w
        
        #=========================================================================================================
        
    def plot(self,color=None,ls='solid',alpha=1,lw=1, n_max=np.inf):
        """
        Plot of the profile apertures.

        Parameters
        ----------
        color : str, optional
            Color of the appeertures in the plot. The default is None.
        ls : str, optional
            Linestyle "solid", "dotted", "dashed", "dashdot", ... The default is 'solid'

        Returns
        -------
        Plot of the apertures

        """
                
        # Paleta de colores
        if color==None:
            from matplotlib.pyplot import cm
            paleta=cm.rainbow(np.linspace(0,1,len(self.ap_r)+1))
        else:
            paleta=[color]*(len(self.ap_r)+1)
        
        # Pintamos las aperturas
        self.ap_c.plot(color=paleta[0],ls=ls, alpha=alpha,lw=lw)
        x, y = self.ap_c.positions
        #plt.text(x,y,0,ha='center',va='bottom')
           
        for i, (ap_r0, ap_l0) in enumerate(zip(self.ap_r, self.ap_l)):
            if i < n_max:
                # Derecha
                ap_r0.plot(color=paleta[i+1],ls=ls, alpha=alpha,lw=lw)
                x, y = ap_r0.positions
                #plt.text(x,y,i+1,ha='center',va='bottom',color=paleta[i+1])
                # Izquierda
                ap_l0.plot(color=paleta[i+1],ls=ls,alpha=alpha,lw=lw)
                x, y = ap_l0.positions
                #plt.text(x,y,i+1,ha='center',va='bottom',color=paleta[i+1])
                
    def saveDS9(self, fname):
        """
        Writes the apertures in a file as DS9 regions

        Parameters
        ----------
        fname : str
            File name.

        Returns
        -------
        A .reg file in the current folder.

        """
        # Cabecera igual para todas las regiones
        l1 = '# Region file format: DS9 version 4.1'
        l2 = 'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
        l3= 'image'
        
        file=open(fname+'.reg', 'w')
        file.write(l1+ "\n")
        file.write(l2+ "\n")
        file.write(l3+ "\n")
        
        # Escribimos la linea para la apertura central
        x,y = self.ap_c.positions
        w = self.ap_c.w
        h = self.ap_c.h
        
        l='box('+str(x)+','+str(y)+','+str(w)+','+str(h)+') # text={0} dash=1'
        file.write(l+ "\n")

        # Ahora escribimos para el resto de las aperturas
        for i, (ap_r0, ap_l0) in enumerate(zip(self.ap_r, self.ap_l)):          
            # Aperturas derechas
            x,y = ap_r0.positions
            w = ap_r0.w
            h = ap_r0.h
        
            l_r='box('+str(x)+','+str(y)+','+str(w)+','+str(h)+') # text={'+str(i+1)+'} dash=1'
            file.write(l_r+ "\n")
            
            # Aperturas izquierdas
            x,y = ap_l0.positions
            w = ap_l0.w
            h = ap_l0.h
        
            l_l='box('+str(x)+','+str(y)+','+str(w)+','+str(h)+') # text={'+str(i+1)+'} dash=1'
            file.write(l_l+ "\n")
    
        file.close()
    
    
    
    def surfbrigth(self, method,subpixels=None,sigma=3):
        """
        Extracts the surface brightness profile from the apertures defined in CentralProfile class.

        Parameters
        ----------
 
        method : interpolation method
            'exact': the exact intersection of the aperture with each pixel is calculated. 
            'center': a pixel is considered to be entirely in or out of the aperture depending on whether its center is in or out of the aperture. 
            'subpixel': pixels are divided into a number of subpixels, which are in or out of the aperture based on their centers. For this method, the number of subpixels needs to be set with the subpixels keyword.
            'center' and 'subpixel', are faster, but with the expense of less precision.
        subpixels : int, optional
            If 'subpixel' method is used.. The default is None.
        sigma : float, optional
            Sigma value for the sigma-clipping rejection algorithm in the statistics. The default is 3.

        Returns
        -------
        astropy.table.Table
            Table with profile values in counts and surface brightness, and their errors.

        """
        
        from astropy.table import Table
        from astropy.stats import sigma_clipped_stats  
        from astropy.stats import sigma_clip
        # Clase auxiliar que define la transformacion de cuantas a brillo superficial         
        c2sb = Counts2Sb(self.zp,self.pixscale,self.A)

        #==================================================================================
        #  El primer elemento tanto a la derecha como a la izquierda es la apertura central
        #==================================================================================
        # =============================================================================
        # Variables espaciales
        # =============================================================================
        x = [0]
        y = [0]
        
        # seccion dentro de la apertura de la posicion central
        sec_data = SecImData(self.ap_c,self.image,method,subpixels=subpixels)
        # media sigma clip
        median = sigma_clipped_stats(sec_data, sigma=sigma)[1]  
        
        # =============================================================================
        # Brillo superficial         
        # =============================================================================
        # Cuentas
        c_r=[median]
        c_l=[median]
        c = [median]
        
        # sum flux in the bin
        sum_ct = np.nansum(sec_data)
        sum_sb = [c2sb.do(sum_ct)]
        
        
        # Magnitudes
        sb_r = [c2sb.do(median)]
        sb_l = [c2sb.do(median)]
        sb = [c2sb.do(median)]
        
        # =============================================================================
        # Errores         
        # =============================================================================        
        # error Poissoniano
        error_poiss=np.sqrt(np.nanmean(np.square(median - sec_data))/np.size(sec_data))      
        # error bin
        error_bin=self.rms_sky/np.sqrt(float(np.size(sec_data)))
        
    
        # Cuentas
        # error total
        err_t = np.sqrt(np.square(error_poiss)+np.square(error_bin))
        c_err = [err_t]
        c_err_r = [err_t]
        c_err_l = [err_t]
 
        # Magnitudes
        # error total
        up = median + err_t
        down = median - err_t
        
        # Pasamos los errores a brillo superficial
        error_down_r=[-2.5*np.log10(median/up)] 
        error_up_r=[-2.5*np.log10(down/median)]

        error_down_l=[-2.5*np.log10(median/up)] 
        error_up_l=[-2.5*np.log10(down/median)]
        
        error_down=[-2.5*np.log10(median/up)] 
        error_up=[-2.5*np.log10(down/median)]
        
        # =============================================================================
        #  Bin size in pixels       
        # =============================================================================
        bin_size = [float(np.size(sec_data))]
        
        #==================================================================================
        # El resto de las aperturas
        #==================================================================================
        for ap_r0, ap_l0 in zip(self.ap_r,self.ap_l):
            # Variables espaciales
            x.append(ap_r0.positions[0]-self.ap_c.positions[0])
            y.append(ap_r0.positions[1]-self.ap_c.positions[1])

            
            # Aperturas derechas
            sec_data_r = SecImData(ap_r0,self.image,method,subpixels=subpixels) 
            # media sigma clip
            median_r= sigma_clipped_stats(sec_data_r, sigma=sigma)[1]


            # Aperturas izquierdas
            sec_data_l = SecImData(ap_l0,self.image,method,subpixels=subpixels)
            # media sigma clip
            median_l= sigma_clipped_stats(sec_data_l, sigma=sigma)[1]
            
            mean_rl = np.nanmean([median_r,median_l])

            # =============================================================================
            # Brillo superficial         
            # =============================================================================
            # Cuentas
            c_r.append(median_r)
            c_l.append(median_l)
            c.append(mean_rl)
            
            # Magnitudes
            sb_r.append(c2sb.do(median_r))
            sb_l.append(c2sb.do(median_l))
            sb.append(c2sb.do(mean_rl))
        
            # sum flux in the bin
            sum_ct_r = np.nansum(sec_data_r)
            sum_sb_r = c2sb.do(sum_ct_r)
            sum_ct_l = np.nansum(sec_data_l)
            sum_sb_l = c2sb.do(sum_ct_l)
            sum_ct = np.nanmean([sum_ct_r, sum_ct_l])
            sum_sb.append(c2sb.do(sum_ct)) 
            
            
            #print(np.nansum(sec_data_r), np.nansum(sigma_clip(sec_data_r, sigma=3)) )
            
            
            # =============================================================================
            # Errores         
            # =============================================================================

            # error Poissoniano
            error_poiss_r = np.sqrt(np.nanmean(np.square(median_r - sec_data_r))/np.size(sec_data_r)) 
            error_poiss_l = np.sqrt(np.nanmean(np.square(median_l - sec_data_l))/np.size(sec_data_l))
            error_poiss = np.sqrt(np.square(error_poiss_r)+np.square(error_poiss_l))
            
            #error bin
            error_bin_r = self.rms_sky/np.sqrt(float(np.size(sec_data_r)))
            error_bin_l = self.rms_sky/np.sqrt(float(np.size(sec_data_l)))
            error_bin = np.sqrt(np.square(error_bin_r)+np.square(error_bin_l))


            # Cuentas
            # error total y a cada lado (derecha/izda)
            err_r = np.sqrt(np.square(error_poiss_r) + np.square(error_bin_r))
            c_err_r.append(err_r)
            
            err_l = np.sqrt(np.square(error_poiss_l) + np.square(error_bin_l))
            c_err_l.append(err_l)
            
            err_t = np.sqrt(np.square(error_poiss) + np.square(error_bin))
            c_err.append(err_t)

            # Brillo superficial
            # error total y a cada lado (derecha/izda)
            up_r = median_r + err_r
            down_r = median_r - err_r
            
            up_l = median_l + err_l
            down_l = median_l - err_l
            
            up = mean_rl + err_t
            down = mean_rl - err_t
            
            # Pasamos los errores a brillo superficial
            error_down_r.append(-2.5*np.log10(median_r/up_r)) 
            error_up_r.append(-2.5*np.log10(down_r/median_r))

            error_down_l.append(-2.5*np.log10(median_l/up_l)) 
            error_up_l.append(-2.5*np.log10(down_l/median_l))

            error_down.append(-2.5*np.log10(mean_rl/up)) 
            error_up.append(-2.5*np.log10(down/mean_rl))
                        
            # =============================================================================
            #    Bin size in pixels         
            # =============================================================================
            bin_size.append(float(np.size(sec_data_r)))     
            
            
        # =============================================================================
        # Reemplazamos los nan por 0 
        # =============================================================================
        error_up = np.nan_to_num(error_up)
        error_down = np.nan_to_num(error_down)

        error_down_r = np.nan_to_num(error_down_r)
        error_up_r = np.nan_to_num(error_up_r)

        error_down_l = np.nan_to_num(error_down_l)
        error_up_l = np.nan_to_num(error_up_l)


        #==================================================================================
        # Definimos el diccionario de la variable espacial y le damos el valor en arcsec
        #==================================================================================
        self.param['horizontal'] = {'spatial' : ('r', x)}
        self.param['vertical'] = {'spatial' : ('z', y)}

        spatial_name = self.param[self.orientation]['spatial'][0]
        spatial = np.array(self.param[self.orientation]['spatial'][1]) * self.pixscale
        
        bins_width = self.w * self.pixscale
        
        vec = [spatial,
               sb, 
               error_up,
               error_down,
               sb_r, 
               sb_l, 
               error_up_r,
               error_down_r,
               error_up_l,
               error_down_l,
               c,
               c_err, 
               c_r,
               c_l,
               c_err_r,
               c_err_l,
               bin_size,
               sum_sb,
               bins_width]        
             
        names = [spatial_name,
                 'sb', 
                 'err bar_up', 
                 'err bar_down', 
                 'sb right',
                 'sb left',
                 'err right bar_up',
                 'err right bar_down',
                 'err left bar_up',
                 'err left bar_down',
                 'ct', 
                 'ct err', 
                 'ct right',
                 'ct left',
                 'ct err right',
                 'ct err left',
                 'bin size pix',
                 'sum flux in sb',
                 'bins width']

        return Table(vec, names=names, meta={'Surface Brightnes': 'table info, u.arcsec'})



#==============================================================================
# Perfiles radiales desplazados        
#==============================================================================
            
class ShiftedProfile:
    
    def __init__(self,image,Xc,Yc,nbins,npix,height,zp,pixscale,delta,rms_sky=0,A=0,orientation='horizontal'):
        
        '''
        Class definition to extract a surface brightness profile with an offset from the main axis of the galaxies using rectangular apertures over an image. 
        Here the apertures are built.    
        
        
        Parameters
        ----------
        image : 2D numpy array
            Data where to extract the profile..
        Xc : float
            X center coordinate in pixels.
        Yc : float
            Y center coordinate in pixels.
        nbins : TYPE
            Number of bins of the profile.
        npix : TYPE
            The full width of the profile in pixels. Number of pixels up to where extract the profile.
        height : float
            The full height of the rectangles in pixels.
        zp : float
            Zero point [mag/arcsec^2].
        delta : float
            Offset in pixels from the galaxy center where to extract the profile.
        rms_sky : float
            RMS of the sky to include as a source of error. 
        A : float
            Galactic absorption or extinction for a given band. 
        orientation : str, optional
            Orientation of the profile.'horizontal' is along the galaxy major axis; 'vertical' if along the minor axis.
            The default is 'horizontal'. 

        Returns
        -------
        None.

        '''
        self.image=image
        self.Xc = Xc
        self.Yc = Yc
        self.nbins = nbins
        self.npix = npix
        self.height = height
        self.zp = zp
        self.A = A
        self.rms_sky = rms_sky
        self.orientation=orientation
        self.delta=delta
        self.pixscale=pixscale
        # Instanciamos a la clase principal
        self.up = CentralProfile(image,Xc,Yc,nbins,npix,height,zp,pixscale,A,rms_sky,orientation,delta)
        self.down = CentralProfile(image,Xc,Yc,nbins,npix,height,zp,pixscale,A,rms_sky,orientation,-delta)

    def plot(self,color=None,ls='solid',alpha=1,lw=1,n_max=np.inf):
        '''
        Plot of the profile apertures.

        Parameters
        ----------
        color : str, optional
            Color of the appeertures in the plot. The default is None.
        ls : str, optional
            Linestyle "solid", "dotted", "dashed", "dashdot", ... The default is 'solid'

        Returns
        -------
        Plot of the apertures
        '''
        
        self.up.plot(color=color,ls=ls,alpha=alpha,lw=lw,n_max=n_max)
        self.down.plot(color=color,ls=ls,alpha=alpha,lw=lw,n_max=n_max)
    
    def saveDS9(self, fname):
        '''
        Writes the apertures in a file as DS9 regions

        Parameters
        ----------
        fname : str
            File name.

        Returns
        -------
        A .reg file in the current folder.

        '''
        self.up.saveDS9(fname+'_up')
        self.down.saveDS9(fname+'_down')        
        

    def surfbrigth(self,method,subpixels=None,sigma=3):
        '''
        Extracts the surface brightness profile from the apertures defined in ShiftedProfile class.

        Parameters
        ----------
 
        method : interpolation method
            'exact': the exact intersection of the aperture with each pixel is calculated. 
            'center': a pixel is considered to be entirely in or out of the aperture depending on whether its center is in or out of the aperture. 
            'subpixel': pixels are divided into a number of subpixels, which are in or out of the aperture based on their centers. For this method, the number of subpixels needs to be set with the subpixels keyword.
            'center' and 'subpixel', are faster, but with the expense of less precision.
        subpixels : int, optional
            If 'subpixel' method is used.. The default is None.
        sigma : float, optional
            Sigma value for the sigma-clipping rejection algorithm in the statistics. The default is 3.

        Returns
        -------
        astropy.table.Table
            Table with profile values in counts and surface brightness, and their errors.

        '''
        
        # Clase auxiliar que define la transformacion de cuentas a brillo superficial         
        c2sb = Counts2Sb(self.zp,self.pixscale,self.A)

        t_1=self.up.surfbrigth(method,subpixels,sigma)
        t_2=self.down.surfbrigth(method,subpixels,sigma)
        
        # Cogemos el Nombre y valor de la coordenada espacial directamente de una de las tablas (son las mismas pada ambas)
        spatial=t_1.columns[0]
        spatial_name=t_1.colnames[0]
       
        # Aperturas arriba/abajo
        # Cuentas
        ct_up=t_1['ct']
        ct_down=t_2['ct']
        # Surface brightness arriba/abajo
        sb_up = t_1['sb']
        sb_down = t_2['sb']
        

        # Aperturas derechas
        # Cuentas
        ct_up_r=t_1['ct right']
        ct_down_r=t_2['ct right']
        ct_r = np.nanmean([ct_up_r,ct_down_r],0)
        # Surface brightness derechas
        sb_r = c2sb.do(ct_r) 
        sb_up_r = t_1['sb right']
        sb_down_r = t_2['sb right']
        
        
        
        # Aperturas izquierdas
        # Cuentas
        ct_up_l = t_1['ct left']
        ct_down_l = t_2['ct left']
        ct_l = np.nanmean([ct_up_l,ct_down_l],0)
        # Surface brightness izquierdas
        sb_l = c2sb.do(ct_l) 
        sb_up_l = t_1['sb left']
        sb_down_l = t_2['sb left']
        
        
        # Media de TODAS las aperturas
        # Cuentas
        ct_mean = np.nanmean([ct_up_r, ct_down_r, ct_up_l, ct_down_l],0)
        # Surface brightness de todas las aperturas
        sb_mean = c2sb.do(ct_mean) 

        
        # =============================================================================
        # Errores
        # =============================================================================
        
        # Error en cuentas
        # Arriba/abajo
        ct_err_up = t_1['ct err']
        ct_err_down = t_2['ct err']        

        # Derechas
        ct_err_up_r = t_1['ct err right']
        ct_err_down_r = t_2['ct err right']  
        ct_err_r = np.sqrt(np.square(ct_err_up_r)+np.square(ct_err_down_r))

        # Izquierdas
        ct_err_up_l = t_1['ct err left']
        ct_err_down_l = t_2['ct err left']   
        ct_err_l = np.sqrt(np.square(ct_err_up_l)+np.square(ct_err_down_l))

        # Error TOTAL en cuentas
        ct_err = np.sqrt(np.square(ct_err_up)+np.square(ct_err_down))        
        
        
        #######################################################################
        # Error brillo superficial
        # Arriba/abajo
        err_up_barDown = - 2.5*np.log10(ct_up/(ct_up + ct_err_up))
        err_up_barUp = - 2.5*np.log10((ct_up - ct_err_up)/ct_up)
        
        err_down_barDown = - 2.5*np.log10(ct_down/(ct_down + ct_err_down))
        err_down_barUp = - 2.5*np.log10((ct_down - ct_err_down)/ct_down)

        # Derechas
        err_up_r_barDown = - 2.5*np.log10(ct_up_r/(ct_up_r + ct_err_up_r))
        err_up_r_barUp = - 2.5*np.log10((ct_up_r - ct_err_up_r)/ct_up_r)

        err_down_r_barDown = - 2.5*np.log10(ct_down_r/(ct_down_r + ct_err_down_r))
        err_down_r_barUp = - 2.5*np.log10((ct_down_r - ct_err_down_r)/ct_down_r)

        err_r_barDown = - 2.5*np.log10(ct_r/(ct_r + ct_err_r))
        err_r_barUp = - 2.5*np.log10((ct_r - ct_err_r)/ct_r)        

        # Izquierdas
        err_up_l_barDown = - 2.5*np.log10(ct_up_l/(ct_up_l + ct_err_up_l))
        err_up_l_barUp = - 2.5*np.log10((ct_up_l - ct_err_up_l)/ct_up_l)

        err_down_l_barDown = - 2.5*np.log10(ct_down_l/(ct_down_l + ct_err_down_l))
        err_down_l_barUp = - 2.5*np.log10((ct_down_l - ct_err_down_l)/ct_down_l)

        err_l_barDown = - 2.5*np.log10(ct_l/(ct_l + ct_err_l))
        err_l_barUp = - 2.5*np.log10((ct_l - ct_err_l)/ct_l)          
        
        # Error TOTAL en brillo superficial
        err_barDown = - 2.5*np.log10(ct_mean/(ct_mean + ct_err))
        err_barUp = - 2.5*np.log10((ct_mean - ct_err)/ct_mean)
        
        

        # =============================================================================
        #  Tabla con los datos  
        # =============================================================================
        from astropy.table import Table
        vec = [spatial,
               sb_mean,
               err_barUp,
               err_barDown,
               sb_up,
               sb_down,
               err_up_barUp,
               err_up_barDown,
               err_down_barUp,
               err_down_barDown,
               sb_r, 
               sb_l, 
               err_r_barUp, 
               err_r_barDown, 
               err_l_barUp, 
               err_l_barDown,
               sb_up_r, 
               sb_down_r, 
               err_up_r_barUp, 
               err_up_r_barDown, 
               err_down_r_barUp, 
               err_down_r_barDown,
               sb_up_l, 
               sb_down_l, 
               err_up_l_barUp, 
               err_up_l_barDown, 
               err_down_l_barUp, 
               err_down_l_barDown,
               ct_mean, 
               ct_err,
               ct_up, 
               ct_down, 
               ct_err_up, 
               ct_err_down,
               ct_r, 
               ct_l, 
               ct_err_r, 
               ct_err_l,
               ct_up_r, 
               ct_down_r, 
               ct_err_up_r, 
               ct_err_down_r,
               ct_up_l, 
               ct_down_l, 
               ct_err_up_l, 
               ct_err_down_l]
                
        names = [spatial_name, 
                 'sb', 
                 'err bar_up', 
                 'err bar_down',
                 'sb up', 
                 'sb down', 
                 'err up bar_up', 
                 'err up bar_down', 
                 'err down bar_up', 
                 'err down bar_down',
                 'sb right', 
                 'sb left', 
                 'err right bar_up', 
                 'err right bar_down', 
                 'err left bar_up', 
                 'err left bar_down',
                 'sb up right', 
                 'sb down right', 
                 'err up right bar_up', 
                 'err up right bar_down', 
                 'err down right bar_up', 
                 'err down right bar_down',
                 'sb up left', 
                 'sb down left', 
                 'err up left bar_up', 
                 'err up left bar_down', 
                 'err down left bar_up', 
                 'err down left bar_down',
                 'ct',
                 'ct err',
                 'ct up', 
                 'ct down', 
                 'ct err up', 
                 'ct err down',
                 'ct right', 
                 'ct left', 
                 'ct err right', 
                 'ct err left',
                 'ct up right', 
                 'ct down right', 
                 'ct err up right', 
                 'ct err down right',
                 'ct up left', 
                 'ct down left', 
                 'ct err up left', 
                 'ct err down left']
       
        return Table(vec, names=names, meta={'Surface Brighnes': 'table info, u.arcsec'})
 
# =============================================================================
# Proyeccion sobre una nueva imagen
# =============================================================================


class ProjectCentralProfile:
    def __init__(self, image, profile, wcs_ori, wcs, zp, pixscale, A=0,rms_sky=0,dx=0, dy=0):

        # Diccionario para todos los parametros segun la apertura sea vertical u horizontal.
        self.param={'horizontal': {'signo' : (-1,1)},            
                     'vertical': {'signo' : (1,-1)}}

        # Geometric parameters
        self.orientation = profile.orientation
        
        # WCS
        self.wcs_ori = wcs_ori 
        self.wcs = wcs
        
        # Data parametes
        self.image=image
        self.zp = zp
        self.A = A
        self.rms_sky = rms_sky
        self.pixscale = pixscale
        self.profile = profile

        #=======================================================================================================
        #  Apertures definition
        #=======================================================================================================
        # Leemos las aperturas de la clase original y las pasamos a las coordenadas en pixel de la otra imagen
        self.ap_c = profile.ap_c.to_sky(wcs_ori).to_pixel(wcs)
        self.ap_c.positions[0] = self.ap_c.positions[0] - dx
        self.ap_c.positions[1] = self.ap_c.positions[1] - dy
        
        self.ap_r=[] # Aperturas hacia la derecha
        self.ap_l=[] # Aperuras hacia la izquierda

        for i in profile.ap_l:
            ap_l = i.to_sky(wcs_ori).to_pixel(wcs)
            ap_l.positions[0] = ap_l.positions[0] - dx
            ap_l.positions[1] = ap_l.positions[1] - dy
            self.ap_l.append(ap_l)            
            
        for i in profile.ap_r:
            ap_r = i.to_sky(wcs_ori).to_pixel(wcs)
            ap_r.positions[0] = ap_r.positions[0] - dx
            ap_r.positions[1] = ap_r.positions[1] - dy
            self.ap_r.append(ap_r)            

         
        
        #=========================================================================================================
        
    def plot(self,color=None,ls='solid'):
        """
        Plot of the profile apertures.

        Parameters
        ----------
        color : str, optional
            Color of the appeertures in the plot. The default is None.
        ls : str, optional
            Linestyle "solid", "dotted", "dashed", "dashdot", ... The default is 'solid'.
        cent_ap : bool, optional
            Plot a mark in the central position of the aperture. The default is False.

        Returns
        -------
        Plot of the apertures

        """
                
        # Paleta de colores
        if color==None:
            from matplotlib.pyplot import cm
            paleta=cm.rainbow(np.linspace(0,1,len(self.ap_r)+1))
        else:
            paleta=[color]*(len(self.ap_r)+1)
        
        # Pintamos las aperturas
        self.ap_c.plot(color=paleta[0],ls=ls)
        x, y = self.ap_c.positions
        plt.text(x,y,0,ha='center',va='bottom')
           
        for i, (ap_r0, ap_l0) in enumerate(zip(self.ap_r, self.ap_l)):
            # Derecha
            ap_r0.plot(color=paleta[i+1],ls=ls)
            x, y = ap_r0.positions
            plt.text(x,y,i+1,ha='center',va='bottom',color=paleta[i+1])
            # Izquierda
            ap_l0.plot(color=paleta[i+1],ls=ls)
            x, y = ap_l0.positions
            plt.text(x,y,i+1,ha='center',va='bottom',color=paleta[i+1])
                
    def saveDS9(self, fname):
        """
        Writes the apertures in a file as DS9 regions

        Parameters
        ----------
        fname : str
            File name.

        Returns
        -------
        A .reg file in the current folder.

        """
        # Cabecera igual para todas las regiones
        l1 = '# Region file format: DS9 version 4.1'
        l2 = 'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
        l3= 'image'
        
        file=open(fname+'.reg', 'w')
        file.write(l1+ "\n")
        file.write(l2+ "\n")
        file.write(l3+ "\n")
        
        # Escribimos la linea para la apertura central
        x,y = self.ap_c.positions
        w = self.ap_c.w
        h = self.ap_c.h
        
        l='box('+str(x)+','+str(y)+','+str(w)+','+str(h)+') # text={0} dash=1'
        file.write(l+ "\n")

        # Ahora escribimos para el resto de las aperturas
        for i, (ap_r0, ap_l0) in enumerate(zip(self.ap_r, self.ap_l)):          
            # Aperturas derechas
            x,y = ap_r0.positions
            w = ap_r0.w
            h = ap_r0.h
        
            l_r='box('+str(x)+','+str(y)+','+str(w)+','+str(h)+') # text={'+str(i+1)+'} dash=1'
            file.write(l_r+ "\n")
            
            # Aperturas izquierdas
            x,y = ap_l0.positions
            w = ap_l0.w
            h = ap_l0.h
        
            l_l='box('+str(x)+','+str(y)+','+str(w)+','+str(h)+') # text={'+str(i+1)+'} dash=1'
            file.write(l_l+ "\n")
    
        file.close()
    
    
    
    def surfbrigth(self, method,subpixels=None,sigma=3):
        """
        Extracts the surface brightness profile from the apertures defined in CentralProfile class.

        Parameters
        ----------
 
        method : interpolation method
            'exact': the exact intersection of the aperture with each pixel is calculated. 
            'center': a pixel is considered to be entirely in or out of the aperture depending on whether its center is in or out of the aperture. 
            'subpixel': pixels are divided into a number of subpixels, which are in or out of the aperture based on their centers. For this method, the number of subpixels needs to be set with the subpixels keyword.
            'center' and 'subpixel', are faster, but with the expense of less precision.
        subpixels : int, optional
            If 'subpixel' method is used.. The default is None.
        sigma : float, optional
            Sigma value for the sigma-clipping rejection algorithm in the statistics. The default is 3.

        Returns
        -------
        astropy.table.Table
            Table with profile values in counts and surface brightness, and their errors.

        """
        
        from astropy.table import Table
        from astropy.stats import sigma_clipped_stats        
        # Clase auxiliar que define la transformacion de cuentas a brillo superficial         
        c2sb = Counts2Sb(self.zp,self.pixscale,self.A)

        #==================================================================================
        #  El primer elemento tanto a la derecha como a la izquierda es la apertura central
        #==================================================================================
        # =============================================================================
        # Variables espaciales
        # =============================================================================
        x = [0]
        y = [0]
        
        # seccion dentro de la apertura de la posicion central
        sec_data = SecImData(self.ap_c,self.image,method,subpixels=subpixels)
        # media sigma clip
        median = sigma_clipped_stats(sec_data, sigma=sigma)[1]  
        
        
        # =============================================================================
        # Brillo superficial         
        # =============================================================================
        # Cuentas
        c_r=[median]
        c_l=[median]
        c = [median]
        
        # Magnitudes
        sb_r = [c2sb.do(median)]
        sb_l = [c2sb.do(median)]
        sb = [c2sb.do(median)]
        
        # =============================================================================
        # Errores         
        # =============================================================================        
        # error Poissoniano
        error_poiss=np.sqrt(np.nanmean(np.square(median - sec_data))/np.size(sec_data))      
        # error bin
        error_bin=self.rms_sky/np.sqrt(float(np.size(sec_data)))
        
    
        # Cuentas
        # error total
        err_t = np.sqrt(np.square(error_poiss)+np.square(error_bin))
        c_err = [err_t]
 
        # Magnitudes
        # error total
        up = median + err_t
        down = median - err_t
        
        # Pasamos a magnitudes
        error_down=[-2.5*np.log10(median/up)] 
        error_up=[-2.5*np.log10(down/median)]
        
        #==================================================================================
        # El resto de las aperturas
        #==================================================================================
        for ap_r0, ap_l0 in zip(self.ap_r,self.ap_l):
            # Variables espaciales
            # Proyectamos las coordenadas espaciales al sistema de coordenadas original.
            x0 = ap_r0.to_sky(self.wcs).to_pixel(self.wcs_ori).positions[0]- self.ap_c.to_sky(self.wcs).to_pixel(self.wcs_ori).positions[0]
            y0 = ap_r0.to_sky(self.wcs).to_pixel(self.wcs_ori).positions[1]- self.ap_c.to_sky(self.wcs).to_pixel(self.wcs_ori).positions[1]
            
            x.append(x0)
            y.append(y0)

            
            # Aperturas derechas
            sec_data_r = SecImData(ap_r0,self.image,method,subpixels=subpixels) 
            # media sigma clip
            median_r= sigma_clipped_stats(sec_data_r, sigma=sigma)[1]

            # Aperturas izquierdas
            sec_data_l = SecImData(ap_l0,self.image,method,subpixels=subpixels)
            # media sigma clip
            median_l= sigma_clipped_stats(sec_data_l, sigma=sigma)[1]
            
            mean_rl = np.nanmean([median_r,median_l])

            # =============================================================================
            # Brillo superficial         
            # =============================================================================
            # Cuentas
            c_r.append(median_r)
            c_l.append(median_l)
            c.append(mean_rl)
            
            # Magnitudes
            sb_r.append(c2sb.do(median_r))
            sb_l.append(c2sb.do(median_l))
            sb.append(c2sb.do(mean_rl))
        
            # =============================================================================
            # Errores         
            # =============================================================================

            # error Poissoniano
            error_poiss_r=np.sqrt(np.nanmean(np.square(median_r - sec_data_r))/np.size(sec_data_r)) 
            error_poiss_l=np.sqrt(np.nanmean(np.square(median_l - sec_data_l))/np.size(sec_data_l))
            error_poiss=np.sqrt(np.square(error_poiss_r)+np.square(error_poiss_l))
            
            #error bin
            error_bin_r=self.rms_sky/np.sqrt(float(np.size(sec_data_r)+np.size(sec_data_r)))
            error_bin_l=self.rms_sky/np.sqrt(float(np.size(sec_data_l)+np.size(sec_data_l)))
            error_bin=np.sqrt(np.square(error_bin_r)+np.square(error_bin_l))


            # Cuentas
            # error total
            err_t = np.sqrt(np.square(error_poiss) + np.square(error_bin))
            c_err.append(err_t)

            # Magnitudes
            # error total
            up = mean_rl + err_t
            down = mean_rl - err_t
            
            # Pasamos los errores a magnitudes
            error_down.append(- 2.5*np.log10(mean_rl/up)) 
            error_up.append(- 2.5*np.log10(down/mean_rl))
            
        # =============================================================================
        # Reemplazamos los nan por 0 
        # =============================================================================
        error_up = np.nan_to_num(error_up)
        error_down = np.nan_to_num(error_down)

        #==================================================================================
        # Definimos el diccionario de la variable espacial y le damos el valor en arcsec
        #==================================================================================
        self.param['horizontal'] = {'spatial' : ('r', x)}
        self.param['vertical'] = {'spatial' : ('z', y)}

        spatial_name = self.param[self.orientation]['spatial'][0]
        spatial = np.array(self.param[self.orientation]['spatial'][1]) * self.profile.pixscale # Aqui el valor de la escala de pixel del sistema original
        
        vec = [spatial, sb_r, sb_l,sb,error_up,error_down,c_r, c_l,c,c_err]
        names = [spatial_name, 'sb right', 'sb left','sb','err up','err down', 'ct right', 'ct left','ct','ct err']
        
        return Table(vec, names=names, meta={'Surface Brightnes': 'table info'})



#==============================================================================
# Perfiles radiales desplazados        
#==============================================================================
            
class ProjectedShiftedProfile:
    
    def __init__(self, image, profile, wcs_ori, wcs, zp, pixscale, A=0, rms_sky=0, dx=0, dy=0):
        
        """
        Class definition to extract a surface brightness profile with an offset from the main axis of the galaxies using rectangular apertures over an image. 
        Here the apertures are built.    
        
        
        Parameters
        ----------
        image : 2D numpy array
            Data where to extract the profile..
        Xc : float
            X center coordinate in pixels.
        Yc : float
            Y center coordinate in pixels.
        nbins : TYPE
            Number of bins of the profile.
        npix : TYPE
            The full width of the profile in pixels. Number of pixels up to where extract the profile.
        height : float
            The full height of the rectangles in pixels.
        zp : float
            Zero point [mag/arcsec^2].
        delta : float
            Offset in pixels from the galaxy center where to extract the profile.
        rms_sky : float
            RMS of the sky to include as a source of error. 
        A : float
            Galactic absorption or extinction for a given band. 
        orientation : str, optional
            Orientation of the profile.'horizontal' is along the galaxy major axis; 'vertical' if along the minor axis.
            The default is 'horizontal'. 

        Returns
        -------
        None.

        """
        
        # Diccionario para todos los parametros segun la apertura sea vertical u horizontal.
        self.param={'horizontal': {'signo' : (-1,1)},            
                     'vertical': {'signo' : (1,-1)}}

        # Geometric parameters
        self.orientation = profile.orientation
        
        # WCS
        self.wcs_ori = wcs_ori 
        self.wcs = wcs
        
        # Data parametes
        self.image=image
        self.zp = zp
        self.A = A
        self.rms_sky = rms_sky
        self.pixscale = pixscale
        self.profile = profile

        #=======================================================================================================
        #  Apertures definition
        #=======================================================================================================
        
        up_ori = profile.up
        down_ori = profile.down
        
        self.up = ProjectCentralProfile(image, up_ori, wcs_ori, wcs, zp, pixscale, A=A,rms_sky=rms_sky,dx=dx, dy=dy)
        self.down = ProjectCentralProfile(image, down_ori, wcs_ori, wcs, zp, pixscale, A=A,rms_sky=rms_sky,dx=dx, dy=dy)

    def plot(self,color=None,ls='solid'):
        '''
        Plot of the profile apertures.

        Parameters
        ----------
        color : str, optional
            Color of the appeertures in the plot. The default is None.
        ls : str, optional
            Linestyle "solid", "dotted", "dashed", "dashdot", ... The default is 'solid'.
        cent_ap : bool, optional
            Plot a mark in the central position of the aperture. The default is False.

        Returns
        -------
        Plot of the apertures
        '''
        
        self.up.plot(color=color,ls=ls)
        self.down.plot(color=color,ls=ls)
    
    def saveDS9(self, fname):
        '''
        Writes the apertures in a file as DS9 regions

        Parameters
        ----------
        fname : str
            File name.

        Returns
        -------
        A .reg file in the current folder.

        '''
        self.up.saveDS9(fname+'_up')
        self.down.saveDS9(fname+'_down')        
        

    def surfbrigth(self,method,subpixels=None,sigma=3):
        '''
        Extracts the surface brightness profile from the apertures defined in ShiftedProfile class.

        Parameters
        ----------
 
        method : interpolation method
            'exact': the exact intersection of the aperture with each pixel is calculated. 
            'center': a pixel is considered to be entirely in or out of the aperture depending on whether its center is in or out of the aperture. 
            'subpixel': pixels are divided into a number of subpixels, which are in or out of the aperture based on their centers. For this method, the number of subpixels needs to be set with the subpixels keyword.
            'center' and 'subpixel', are faster, but with the expense of less precision.
        subpixels : int, optional
            If 'subpixel' method is used.. The default is None.
        sigma : float, optional
            Sigma value for the sigma-clipping rejection algorithm in the statistics. The default is 3.

        Returns
        -------
        astropy.table.Table
            Table with profile values in counts and surface brightness, and their errors.

        '''
        
        # Clase audiliar que define la transformacion de cuantas a brillo superficial         
        c2sb = Counts2Sb(self.zp,self.pixscale,self.A)

        t_1=self.up.surfbrigth(method,subpixels,sigma)
        t_2=self.down.surfbrigth(method,subpixels,sigma)
        
        # Cogemos el Nombre y valor de la coordenada espacial directamente de una de las tablas (son las mismas pada ambas)
        spatial=t_1.columns[0]
        spatial_name=t_1.colnames[0]

        # Aperturas derechas
        # Cuentas
        c1_r=t_1['ct right']
        c2_r=t_2['ct right']
        c_r = np.nanmean([c1_r,c2_r],0)
        
        # Magnitudes
        sb_r = c2sb.do(c_r) 
        
        # Aperturas izquierdas
        # Cuentas
        c1_l = t_1['ct left']
        c2_l = t_2['ct left']
        c_l = np.nanmean([c1_l,c2_l],0)

        # Magnitudes 
        sb_l = c2sb.do(c_l) 
        
        
        # Media de todas las aperturas
        # Cuentas
        c_mean = np.nanmean([c1_r,c2_r,c1_l,c2_l],0)
        # Magnitudes
        sb_mean = c2sb.do(c_mean) 
        
        # =============================================================================
        # Errores
        # =============================================================================
        c1_error = t_1['ct err']
        c2_error = t_2['ct err']        
        
        # Suma cuadratica
        err_t = np.sqrt(np.square(c1_error)+np.square(c2_error))
        
        # Error cuentas
        c_err = np.copy(err_t)

        up = c_mean + err_t
        down = c_mean - err_t
        
        error_down = - 2.5*np.log10(c_mean/up)
        error_up = - 2.5*np.log10(down/c_mean)

        from astropy.table import Table
        vec = [spatial, sb_r, sb_l,sb_mean,error_up,error_down, c_r,c_l,c_mean,c_err]
        
        
        names = [spatial_name, 'sb right', 'sb left','sb','err up','err down','ct right', 'ct left','ct','ct err']

        return Table(vec, names=names, meta={'Surface Brighnes': 'table info'})



    
    
       

       
#********************************************************************************
# ELLIPTICAL RADIAL PROFILE 
#********************************************************************************

class EllipticalProfile:

    def __init__(self, data, Xc, Yc, sma, eps, pa, astep=0.1, linear_growth=False, fix_center=False, fix_pa=False, fix_eps=False, ell_threshold=0.1):
        '''
        Class definition to extract elliptical isophotes values over an image and obtain a profile.
        Here it is defined the geometry of the initial ellipse

        Parameters
        ----------
        data : 2D numpy array
            Image with the data where to extract the profile.
        Xc : float
            X center coordinate in pixels.
        Yc : float
            Y center coordinate in pixels.
        sma : float
            Semi-major axis of the initial ellipse in pixels (but not the smallest).
        eps : float
            The ellipticity of the initial ellipse.
        pa : float
            The position angle (in radians) of the semimajor axis in relation to the postive x axis of the image array (rotating towards the positive y axis). 
            Position angles are defined in the range :math:`0 < PA <= \pi`. Avoid using as starting position angle of 0., since the fit algorithm may not work properly. 
            When the ellipses are such that position angles are near either extreme of the range, noise can make the solution jump back and forth between successive isophotes, by amounts close to 180 degrees.
        astep : float, optional
            The step value for growing/shrinking the semimajor axis. 
            It can be expressed either in pixels (when ``linear_growth=True``) or as a relative value (when ``linear_growth=False``). 
            The default is 0.1.
        linear_growth : bool, optional
            The semimajor axis growing/shrinking mode. The default is 'False'.
        fix_center : bool, optional
            Keep center of ellipse fixed during fit? The default is False.
        fix_pa : bool, optional
            Keep position angle of semi-major axis of ellipse fixed during fit? The default is False.
        fix_eps : bool, optional
            Keep ellipticity of ellipse fixed during fit? The default is False.
        ell_threshold : float, optional
            The threshold for the object centerer algorithm. 
            By lowering this value the object centerer becomes less strict, in the sense that it will accept lower signal-to-noise data. 
            If set to a very large value, the centerer is effectively shut off. In this case, either the geometry information supplied by the ``geometry`` parameter is used as is, or the fit algorithm will terminate prematurely. Note that once the object centerer runs successfully, the (x, y) coordinates in the ``geometry`` attribute (an `~photutils.isophote.EllipseGeometry` instance) are modified in place. 
            The default is 0.1.

        Returns
        -------
        None.

        '''

        self.data = data
        self.Xc = Xc
        self.Yc = Yc
        
        # =============================================================================
        # Instanciamos EllipseGeometry a una variable global de la clase         
        # =============================================================================
        from photutils.isophote import EllipseGeometry     
        self.geometry = EllipseGeometry(self.Xc, self.Yc,sma = sma, eps= eps, pa=pa, astep=astep, linear_growth=linear_growth, fix_center=fix_center, fix_pa=fix_pa, fix_eps=fix_eps)
        
        # =============================================================================
        # Instanciamos Ellipse         
        # =============================================================================
        from photutils.isophote import Ellipse
        self.ellipse = Ellipse(self.data, self.geometry, threshold=ell_threshold)
        
    def find_center(self):
        '''
        Find the center coordinates of the galaxy.

        Returns
        -------
        If the algorithm is successful the (x, y) coordinates in this ~photutils.isophote.EllipseGeometry (i.e., the x0 and y0 attributes) instance will be modified.

        '''
        self.geometry.find_center(self.data)

    
    def fit_image(self, sma0=None, minsma=0.0, maxsma=None, step=0.1, conver=0.05, minit=10, maxit=50, fflag=0.5, maxgerr=0.5, sclip=3.0, nclip=3, integrmode='median', linear=None, maxrit=None, fix_center=False, fix_pa=False, fix_eps=False):
        '''
        Fit multiple isophotes to the image array plus an extra fit for the central region. 
        This method loops over each value of the semimajor axis (sma) length, fitting a single isophote at each sma. 
        The entire set of isophotes is returned in an `~photutils.isophote.IsophoteList` instance.
        More info in: photutils/isophote/ellipse.py package fit_image 

        Parameters
        ----------
        sma0 : float, optional
            The starting value for the semimajor axis length (pixels).
            This value must not be the minimum or maximum semimajor axis length, but something in between. 
            The ``sma0`` value should be selected such that the corresponding isophote has a good signal-to-noise ratio and a clearly defined geometry. 
            If set to `None` (the default), one of two actions will be taken:  
            If a`~photutils.isophote.EllipseGeometry` instance was input to the `~photutils.isophote.Ellipse` constructor, its ``sma`` value will be used.
            Otherwise, a default value of 10. will be used.
        minsma : float, optional
            The minimum value for the semimajor axis length (pixels). 
            The default is 0.0.
        maxsma : float, optional
            The maximum value for the semimajor axis length (pixels).
            When set to `None` (default), the algorithm will increase the semimajor axis until one of several conditions will cause it to stop and revert to fit ellipses with sma < ``sma0``.
        step : float, optional
            The step value used to grow/shrink the semimajor axis length (pixels if ``linear=True``, or a relative value if ``linear=False``). 
            See the ``linear`` parameter. 
            The default is 0.1.
        conver : float, optional
            DESCRIPTION. The default is 0.05.
        minit : int, optional
            DESCRIPTION. The default is 10.
        maxit : int, optional
            DESCRIPTION. The default is 50.
        fflag : TYPE, optional
            DESCRIPTION. The default is 0.7.
        maxgerr : float, optional
            DESCRIPTION. The default is 0.5.
        sclip : float, optional
            DESCRIPTION. The default is 3.0.
        nclip : int, optional
            DESCRIPTION. The default is 0.
        integrmode : {'bilinear', 'nearest_neighbor', 'mean', 'median'}, optional
            DESCRIPTION. The default is 'median'.
        linear : bool, optional
            DESCRIPTION. The default is None.
        maxrit : float, optional
            DESCRIPTION. The default is None.
        fix_center : bool, optional
            DESCRIPTION. The default is False.
        fix_pa : bool, optional
            DESCRIPTION. The default is False.
        fix_eps : bool, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        A `~photutils.isophote.IsophoteList` instance
        A list-like object of `~photutils.isophote.Isophote` instances, sorted by increasing semimajor axis length.

        '''
        self.maxsma = maxsma
        self.isolist = self.ellipse.fit_image(sma0, minsma, self.maxsma, step, conver, minit, maxit, fflag, maxgerr, sclip, nclip, integrmode, linear, maxrit, fix_center, fix_pa, fix_eps)
        
        # =============================================================================
        # Elipse central
        # =============================================================================
        from photutils.isophote.sample import CentralEllipseSample
        from photutils.isophote.fitter import CentralEllipseFitter
        from photutils.isophote import EllipseGeometry

        # central peak position.
        g = EllipseGeometry(self.Xc, self.Yc, 0., 0., 0.)

        # then we build a CentralEllipseSample instance and call the CentralEllipseFitter on it.
        sample = CentralEllipseSample(self.data, 0., geometry=g)
        fitter = CentralEllipseFitter(sample)
        center = fitter.fit()
        # Replace the new fit in the central region
        self.isolist[0] = center
        self.isolist.sort()

    def outer_regions(self,sma_cut,step):
        '''
        Extra fitting of the outer parts of the object if we want to reach very large ragial distances and/or sky regions

        Parameters
        ----------
        sma_cut : float
            Radial position in pix from which we will initiate the fitting for the outer parts
        step : float
            The step value used to grow/shrink the semimajor axis length.
            Use a large value than in 'fit_image' as the SNR here is lower.

        Returns
        -------
        A `~photutils.isophote.IsophoteList` instance
        A list-like object of `~photutils.isophote.Isophote` instances, sorted by increasing semimajor axis length.
        Includes the 'fit_image' fitting + the extended from 'outer_regions'

        '''
        from photutils.isophote import EllipseGeometry
        from photutils.isophote import Ellipse
        import numpy as np

        # Index if the isophote from which we will initiate the fitting for the outer parts
        ind = np.where(self.isolist.sma> sma_cut)[0].min()

        # Initial params from the last ellipse in the previous fitting
        x0 = self.isolist.x0[ind]
        y0 = self.isolist.y0[ind]
        sma = self.isolist.sma[ind]
        eps = self.isolist.eps[ind]
        pa = self.isolist.pa[ind]

        g = EllipseGeometry(x0, y0, sma, eps, pa)
       
        ellipse = Ellipse(self.data, geometry=g)

        outer_isolist = ellipse.fit_image(integrmode='median', step=step, sma0=sma,minsma=sma, fflag=0.2, sclip=3.0, nclip=3, fix_center=True)

        # Add external part of the ellipse fittng
        self.isolist = self.isolist[:ind] + outer_isolist
        
    def run_over(self,new_data, zp, pixscale, rms_sky=0, A=0, dimming=0, k_corr=0):
        '''
        Runs a previous ellipse fitting over a different image using the same ellipses obtained in fit_images (and outer_regions if used)

        Returns
        -------
        Astropy Tablewith surface brightnes values and errors

        '''        
        from photutils.isophote import EllipseSample, Isophote, IsophoteList
        # from photutils.isophote.sample import CentralEllipseSample
        # from photutils.isophote.fitter import CentralEllipseFitter
        # from photutils.isophote import EllipseGeometry

        # =============================================================================
        # Loop over all the isophotes
        # =============================================================================
        new_isolist = []
        
        for iso in self.isolist[1:]:
            # Get the EllipseGeometry of each fitted isophote as determined by 
            # fitting on the high-S/N inage.
            get = iso.sample.geometry
        
            # Sample the new image at the same geometry. Use the same integration 
            # mode so as to ensure the same regions in each image are sampled. 
            # Use the same clipping parameters as used in the high-S/N image.
            sample = EllipseSample(new_data, get.sma, geometry=get,integrmode='median', sclip=3, nclip=3)
            sample.update(fixed_parameters = [False, False, False, False])
        
            # Create an Isophote instance with the sample, and store it in 
            # temporary list. Here we are using '0' as the number of iterations,
            # 'True' for the validity status, and '0' for the stop code. These
            # are in fact arbitrary in this context; you could use anything you
            # like.
            iso_ = Isophote(sample, 0, True, 0)
            new_isolist.append(iso_)
        
        # Build the IsophoteList instance with the result.
        new_isolist = IsophoteList(new_isolist)
        
        # =============================================================================
        # Elipse central
        # =============================================================================
        from photutils.isophote.sample import CentralEllipseSample
        from photutils.isophote.fitter import CentralEllipseFitter
        from photutils.isophote import EllipseGeometry

        # central peak position.
        g = EllipseGeometry(self.Xc, self.Yc, 0., 0., 0.)

        # then we build a CentralEllipseSample instance and call the 
        # CentralEllipseFitter on it.
        sample = CentralEllipseSample(new_data, 0., geometry=g)
        fitter = CentralEllipseFitter(sample)
        center = fitter.fit()
        new_isolist.append(center)
        new_isolist.sort()
        
        return  self.surfbrigth(zp, pixscale, rms_sky=rms_sky, A=A, dimming=dimming, k_corr=k_corr, isolist = new_isolist)

    
    def surfbrigth(self, zp, pixscale, rms_sky=0, A=0, dimming=0, k_corr=0, isolist = None):
        import numpy as np
        from astropy.table import Table
        # If not isolist, take it from the original fit
        if isolist == None:
            isolist = self.isolist
        # Auxiliar clase to convert counts to surface brightness        
        c2sb = Counts2Sb( zp, pixscale, A, dimming, k_corr)
        
        # Read counts from the fit
        sma = isolist.sma
        counts = isolist.intens
        counts_err = np.sqrt(np.square(isolist.int_err/np.sqrt(isolist.npix_e)) + np.square(rms_sky/np.sqrt(isolist.npix_e)))
        
        # SNR
        SNR = counts/isolist.int_err
        
        # Fit parameters
        xc = isolist.x0
        xc_err = isolist.x0_err
        
        yc = isolist.y0
        yc_err = isolist.y0_err
        
        ell = isolist.eps
        ell_err = isolist.ellip_err        
        
        pa = isolist.pa
        pa_err = isolist.pa_err

        b4 = isolist.b4
        b4_err = isolist.b4_err
     
        # Surface brightness value with corrections
        sb = c2sb.do(counts)


        # Errors in magnitudes
        up = counts + counts_err
        down = counts - counts_err

        error_down = -2.5*(np.log10(counts/up))
        error_up = -2.5*(np.log10(down/counts))
        

        
        # =============================================================================
        # Replace nan by 0 
        # =============================================================================
        error_up = np.nan_to_num(error_up)
        error_down = np.nan_to_num(error_down)

        
        # =============================================================================
        # Save a table         
        # =============================================================================

    
        names=['sb','err up','err down','ct','ct err','SNR','sma','xc', 'xc err', 'yc', 'yc err', 'ell', 'ell err', 'pa', 'pa err', 'b4', 'b4 err']
        vec = [sb, error_up, error_down, counts,counts_err,SNR,sma,xc, xc_err, yc, yc_err, ell, ell_err, pa, pa_err, b4, b4_err] 
    
        return Table(vec,names=names, meta={'Surface Brightess': 'table info'})
    

        
    def plot(self,color=None,ls='solid',smamax=np.inf):
        from photutils import EllipticalAperture
        # =============================================================================
        # Read the fir apertures
        # =============================================================================
        print(len(self.isolist))
        x0 = self.isolist.x0
        y0 = self.isolist.y0
        sma = self.isolist.sma
        smna = sma*(1 - self.isolist.eps)
        pa = self.isolist.pa
        # =============================================================================
        # Colors palette/map
        # =============================================================================
        if color==None:
            from matplotlib.pyplot import cm
            paleta=cm.rainbow(np.linspace(0,1,len(x0-1)))
        else:
            paleta=[color]*(len(x0-1))

        for i in range(1,len(x0)):
            if sma[i] <= smamax:
                ap = EllipticalAperture((x0[i], y0[i]), sma[i], smna[i], pa[i])
                ap.plot(color=paleta[i],ls=ls)






# =============================================================================
# Funcion para corregir de NII las imagenes en Halpha
# =============================================================================


def NIIcorr(flux_Ha, M_B):
    '''
    Function to apply the NII correction to Halpha images
    See eq. B1 in Kennicutt+08


    Parameters
    ----------
    flux_Ha : np.array
        Array with Halpha + NII flux values already corrected by galactic absorption. [erg s-1 cm-2] 
    M_B : float
        Absolute magnitude in B band.

    Returns
    -------
    astropy.table.table.Table updated with the Halpha surface brightness values corrected by galactic absorption and NII emission lines contribution.

    '''
    
    if M_B<=-21.:
        nIIcorr = 0.54   # log(NII/Halpha) = 0.54
    else:
        nIIcorr = -0.173*M_B - 3.903
    
    flux_corr = flux_Ha/(1+10**nIIcorr)
    
    return flux_corr

# =============================================================================
# Funcion para pasar magnitudes AB a flujo en uds cgs
# =============================================================================

def magab2fluxcgs(mag_ab, lambda_angs, delta_lambda):
    
    # https://en.wikipedia.org/wiki/AB_magnitude 
    
    flux_jy_nu = 10**(-0.4 * (mag_ab - 8.9))
    flux_cgs_lamda = flux_jy_nu / (3.34e4*(lambda_angs**2))
    flux_cgs = flux_cgs_lamda * delta_lambda
    
    return flux_cgs


# =============================================================================
# Funcion para calcular la SFR de las imagenes en Halpha
# =============================================================================


def SFR_Halpa(flux_cgs, flux_err, DMpc):
    '''
    Funtion to obtain the SFR from Halpha data
    See Kennicutt+2009 and Sect. 5.7 in Sanchez-Gallego+2012:
        
        SFR(M_sun/yr) = 5.5 × 10^−42 L(Hα) (erg s-1) 

    Parameters
    ----------
    flux_cgs : astropy.table.column.Column
        Column of the Halpha flux values in cgs units corrected by galactic absorption and NII emission lines contribution.
    flux_err: astropy.table.column.Column
        Column of the Halpha flux error values in cgs. 
    DMpc : float
        Distance in Mpc to the object

    Returns
    -------
    Array of SFR values in [M_sun/yr].

    '''
    
    Flux = flux_cgs
    
    # Flux in luminosity units: L(Hα) [erg/s]
        # L(Hα)(erg s-1 ) = 4πD^2 * (3.086×10^24)^2 * FHα
    LHa = 4*np.pi*(DMpc**2)*((3.086e24)**2)*Flux
    
    # SFR using the eq. above 
    SFR = LHa * 5.5e-42   #  7.9e-42*LHa [erg s-1] kennicutt 1998
        
    # SFR error is bassically => cte * FHα_err
    SFR_err = (5.5e-42)*4*np.pi*(DMpc**2)*((3.086e24)**2)*flux_err


    vec = [SFR, SFR_err]
    names = ['SFR','SFR err']
        
    sfr_table = Table(vec, names=names, meta={'SFR': 'table info'})

    
    # from astropy.table import hstack
    # sfr_table = hstack([SFR, SFR_err])
    # sfr_table['sb'].name = 'SFR'
    # sfr_table['ct err'].name = 'SFR err'
    
    return sfr_table

            

# =============================================================================
# Funcion para borrar secciones de las mascaras
# =============================================================================


def DeleteMaskRegion(mask, image, vmin= 0, vmax=0.3):
    """ [Interactive] Fuction to delete mask regions by clicking over the regions to take out
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from tkinter import messagebox #tkinter.TkVersion 8.6 
    from astropy.visualization import LogStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    
    plt.ion(),plt.show()

    # Set the limits to ensure that we always plot the same colors in the same masked regions
    vmax_m=np.max(mask)

    remove=True
    plt.figure(1,figsize=(7,7))
    
    while remove==True:
        happy=False
        while happy == False:
            plt.cla()
            plt.imshow(mask,origin='lower',cmap='tab20c', aspect='equal', 
                       vmin=1, vmax=vmax_m, interpolation='none')
            plt.tight_layout()
            plt.show()
        
            mask_new=np.copy(mask)
            image_new=np.copy(image)
            
            # Take points iteratively 
            points=plt.ginput(n=-1,timeout=0,
                              show_clicks=True, mouse_pop=2, mouse_stop=3)
            for point in points:
                coords = np.int_(point)
                value=mask[coords[1]][coords[0]]
                mask_new[mask_new==value]=0
                image_new[mask_new != 0]=np.nan

            
            plt.cla()
            #plt.imshow(mask_new,origin='lower',aspect='auto',vmin=vmin,vmax=vmax)
            
            # Create normalized object
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch())
            
            plt.imshow(image_new,origin='lower',aspect='equal', norm=norm)

            plt.tight_layout()
            plt.show()
            
            happy=messagebox.askyesno("Masks"," Are you happpy with the result?")
            if happy == True:
                mask = np.copy(mask_new)
                remove=messagebox.askyesno("Masks"," Do you want to remove more masked regions?",default='no')
                
    plt.close(1)

    return mask_new  



# =============================================================================
# Funcion para enmascarar una seccion circular en una posicion dada con el cursor
# =============================================================================

def CircularHandMask(image, mask_value = True, vmin= 0, vmax=0.2):
    """ [Interactive] Fuction to add circular mask regions by putting the number of new masks and their radius in pix.
    Then click on the image AS MANY TIMES AS REQUIRED to tell the code where to locate those masks.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from tkinter import messagebox #tkinter.TkVersion 8.6 
    from astropy.visualization import LogStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    from photutils import CircularAperture
    plt.ion()
    shape = np.shape(image)
    
    remove=True
    plt.figure(1,figsize=(7,7))
    
    while remove==True:
        happy=False
        while happy == False:
            plt.cla()
            
            # Create normalizer object
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch())
            plt.imshow(image,origin='lower',aspect='equal', norm=norm)            
            plt.tight_layout()
            plt.show()
        
            image_new=np.copy(image)
            
            # Ask the mask radii (default mask with R=10 pix)
                
            radii=input('masks radii? (default 30 pix): ')
            if radii.isnumeric() == True:
                radii=int(radii)
            else:
                radii=30
            
            # Look for the coords interactively
            coords = plt.ginput(n=-1,timeout=0,show_clicks=True, 
                                mouse_pop=2, mouse_stop=3)
            
            mask = 0
            for coord in coords:
                # Creat a circular aperture at each position
                ap = CircularAperture(coord, radii)
                mask = mask + ap.to_mask(method='center').to_image(shape)
                
            image_new[mask != 0]= mask_value
            
            plt.cla()            
            # Create normalize object
            norm = ImageNormalize(vmin=0., vmax=0.3, stretch=LogStretch())
            
            plt.imshow(image_new,origin='lower',aspect='equal', norm=norm)

            plt.tight_layout()
            plt.show()
            
            happy=messagebox.askyesno("Masks"," Are you happpy with the result?")
            if happy == True:
                image = np.copy(image_new)
                remove=messagebox.askyesno("Masks"," Do you want to continue masking this region?")
                
    plt.close(1)

    return image  
            


def CircularHandMask_new(image, mask_value = True, vmin= 0, vmax=0.2):
    """ [Interactive] Fuction to add circular mask regions by putting the number of new masks and their radius in pix.
    Then click on the image AS MANY TIMES AS REQUIRED to tell the code where to locate those masks.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from tkinter import messagebox #tkinter.TkVersion 8.6 
    from astropy.visualization import LogStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    from photutils import CircularAperture

    shape = np.shape(image)
    
    remove=True
    plt.figure(1,figsize=(7,7))
    
    while remove==True:
        happy=False
        while happy == False:
            plt.cla()
            
            # Create normalizer object
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch())
            plt.imshow(image,origin='lower',aspect='equal', norm=norm)            
            plt.tight_layout()
            plt.show()
        
            image_new = np.copy(image)
            im_zeros = np.zeros_like(image)
            
            # Ask the mask radii (default mask with R=10 pix)
                
            radii=input('masks radii? (default 10 pix): ')
            if radii.isnumeric() == True:
                radii=int(radii)
            else:
                radii=10
            
            # Look for the coords interactively
            coords = plt.ginput(n=-1,timeout=0,show_clicks=True, 
                                mouse_pop=2, mouse_stop=3)
            
            mask = 0
            for coord in coords:
                # Creat a circular aperture at each position
                ap = CircularAperture(coord, radii)
                mask = mask + ap.to_mask(method='center').to_image(shape)
                
            image_new[mask != 0]= mask_value
            im_zeros[mask != 0]= mask_value
            
            plt.cla()            
            # Create normalize object
            norm = ImageNormalize(vmin=0., vmax=0.3, stretch=LogStretch())
            
            plt.close('all')
            plt.imshow(image_new,origin='lower',aspect='equal', norm=norm)

            plt.tight_layout()
            plt.show()
            
            happy=messagebox.askyesno("Masks"," Are you happpy with the result?")
            if happy == True:
                remove=messagebox.askyesno("Masks"," Do you want to continue masking this region?")
                
    plt.close(1)

    return im_zeros  
            

# =============================================================================
# Elliptical masks
# =============================================================================

def EllipticalHandMask(image, mask_value = True, vmin= 0, vmax=0.2):
    """ [Interactive] Function to add elliptical masks 'by hand' just puting how many masks with the same axis values
        Then click on the image AS MANY TIMES AS REQUIRED to tell the code where to locate those masks.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from tkinter import messagebox #tkinter.TkVersion 8.6 
    from astropy.visualization import LogStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    from photutils import EllipticalAperture


    remove=True
    plt.figure(1,figsize=(7,7))
    
    while remove==True:
        happy=False
        while happy == False:
            plt.cla()
            
            # Create normalizer object
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch())
            plt.imshow(image,origin='lower',aspect='equal', norm=norm)            
            plt.tight_layout()
            plt.show()
        
            image_new=np.copy(image)
            
            # Ask how many masks we want to add with the same radii (default 1 mask with a=10 and b=5 pix)
            
            n=input('How many maks? (default 1): ')
            if n.isnumeric() == True:
                n=int(n)
            else:
                n=1
                
            a=input('semi-major axis ? (default 10): ')
            if a.isnumeric() == True:
                a=int(a)
            else:
                a=10

            b=input('semi-minor axis ? (default 5): ')
            if b.isnumeric() == True:
                b=int(b)
            else:
                b=5

            # look for the coords
            coords = plt.ginput(n=n,timeout=0,show_clicks=True)
            
            # Create an elliptical aperture at each position
            ap = EllipticalAperture(coords, a=a, b=b)
            
            shape = np.shape(image)
            
            # Loop over each of the selected positions
            mask = 0
            for i in range(n):
                mask = mask + ap[i].to_mask(method='center').to_image(shape)
            
            image_new[mask != 0]= mask_value
            
            plt.cla()            
            # Create normalizer object
            norm = ImageNormalize(vmin=0., vmax=0.3, stretch=LogStretch())
            
            plt.imshow(image_new,origin='lower',aspect='equal', norm=norm)

            plt.tight_layout()
            plt.show()
            
            happy=messagebox.askyesno("Masks"," Are you happpy with the result?")
            if happy == True:
                image = np.copy(image_new)
                remove=messagebox.askyesno("Masks"," Do you want to continuing masking this region?")
                
    plt.close(1)

    return image  





# =============================================================================
# Seccionar imagen en parches
# =============================================================================

def Patchs(image, size):
    """ Select a patch of an image (data array) of size=size in pixels [square region]
    """
    y_l,x_l = np.shape(image)
    # Create section in the image
    patchs = []
    
    x = np.append(np.arange(0,x_l,size),x_l)
    y = np.append(np.arange(0,y_l,size),y_l)
    
    for i in range(len(x)-1):
        for j in range(len(y)-1):
            sec=slice(y[j],y[j+1]),slice(x[i],x[i+1])  
            patchs.append(sec)
    
    return patchs



# =============================================================================
# Truncation fitting
# =============================================================================

def TruncationFitting(x,y,c=1, a1=1, a2=-1, b1=0,sigma=3,c_min=None, c_max= None):
    from astropy.modeling.models import custom_model
    from astropy.modeling import fitting
    from astropy.stats import sigma_clip


    @custom_model
    def TwoLines(x,c=c, a1=a1, a2=a2, b1=b1):
        import numpy as np
        b2=(a1 - a2)*c+b1
        return(np.concatenate([a1*x[x<=c]+b1,a2*x[x>c]+b2]))

    init = TwoLines()
    init.c.min = c_min
    init.c.max = c_max
    
    fitter = fitting.SLSQPLSQFitter()
    fit = fitting.FittingWithOutlierRemoval(fitter, sigma_clip, niter=3, sigma=sigma)

    return fit(init,x, y)


def TwoLines(x,c=1, a1=1, a2=1, b1=1):
    import numpy as np
    """
    Modelo de dos rectas que se cortan en c, de pendiente a y ordenada en el origen b
    """
    b2=(a1 - a2)*c+b1
    return(np.concatenate([a1*x[x<=c]+b1,a2*x[x>c]+b2]))


# =============================================================================
# LogNorm for imshow in log scale
# =============================================================================

def LogNorm(vmin, vmax):    
    from astropy.visualization import LogStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    
    return ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch())


# =============================================================================
# Funciones para pasar de kpc a arcsec y viceversa
# =============================================================================

# Kpc to Arcsec **************************************************************
def Kpc2Arcsec(x_kpc, DMpc):
    """
    x_Kpc is a value to measure within the image in kpc
    DMpc is the distance to the target in Mpc
    """
    import numpy as np
    
    return np.arctan(x_kpc/(DMpc*1e3))*(180./np.pi)*60*60.


# Arcsec to kpc **************************************************************
def Arcsec2Kpc(x_arcsec, DMpc):
    """
    x_arcsec is a value to measure within the image in arcsec
    DMpc is the distance to the target in Mpc
    """

    import numpy as np
    
    return np.tan((x_arcsec/(60.*60.))*(np.pi/180.))*DMpc*1e3

          
# Pix to kpc **************************************************************
def Pix2Kpc(pix, scale, DMpc=None, z=None, H0=70, Om0=0.3):
    """
    pix is a value to measure within the image in pix, scale is the pix scale in arcsec/pix
    DMpc is the distance to the target in Mpc
    H0 and Om0 are the cosmological parameters
    """
    try:
        """
        Pix to kpc at a given distance in Mpc 
        Angular distance
        """
        import numpy as np
        x_arcsec = pix * scale
        return np.tan((x_arcsec/(60.*60.))*(np.pi/180.))*DMpc*1e3
    except:
        '''
        Pix to kpc at a given redshift AngAperture[kpc] = D_A(z)[Mpc] × R[arcsec]
        Standard cosmological model: H0 = 70 km s−1 Mpc−1, Ωm = 0.3 and ΩΛ = 0.7
        '''
        from astropy.cosmology import FlatLambdaCDM
        from astropy import units as u
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
        d_A = cosmo.angular_diameter_distance(z=z)
        r_arcsec = (pix*scale)*u.arcsec
        return (r_arcsec * d_A).to(u.kpc, u.dimensionless_angles()).value 

# kpc to pix **************************************************************

def Kpc2Pix(kpc, scale, DMpc=None, z=None, H0=70, Om0=0.3):
    try:
        import numpy as np
        x_arcsec = (np.arctan(kpc/(DMpc*1e3))/(np.pi/180.))*60*60
        pix = x_arcsec/scale
        return pix
    except:
        from astropy.cosmology import FlatLambdaCDM
        from astropy import units as u
        
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
        d_A = cosmo.arcsec_per_kpc_proper(z)
        
        x_arcsec = d_A * kpc*u.kpc 
        return (x_arcsec/scale).value




# =============================================================================
#  Desviacion estandar en una region            
# =============================================================================


class RectStat():
    def __init__(self, data, w, h, n=0, Xc=0, Yc=0, theta = 0):
        """ Description
        """
        from photutils import RectangularAperture
        
        if n==0:
            # Valores geometricos de la apertura
            self.ap = [RectangularAperture([Xc,Yc], w=w, h=h, theta=0)]  
        
        else:
            Xc = np.random.rand(n)*(data.shape[1]- 2*w) + w
            Yc = np.random.rand(n)*(data.shape[0]- 2*h) + h
    
            self.ap = []
            
            for i in range(n): 
                # Valores geometricos de la apertura
                self.ap.append(RectangularAperture([Xc[i],Yc[i]], w=w, h=h, theta=0))  
    
            
        self.counts = np.array([])
        
        for ap in self.ap:
            # Seccion de los datos sobre la que haremos la estadistica
            mask=ap.to_mask(method='center') # Geometria de la mascara
            mask_data=mask.data  # Mascara de ceros y unos (segun method) con esa geometria (self.ap)
            sec=mask.cutout(data) # Seccion de los datos de igual forma y tamano que sec 
            sec_weight=sec*mask_data  # Pesamos los datos en funcion del method para determinar la mascara
            self.counts = np.concatenate([self.counts, sec_weight[mask_data>0]])  # Me quedo con la region con los pixeles no enmascarados y ya pesados


    def stat(self, sigma_clip=3):
 
        from astropy.stats import sigma_clipped_stats   
        # mean, median, stddev
        mean, median, stddev = sigma_clipped_stats(self.counts, sigma=sigma_clip, maxiters = 10, mask_value=np.nan)
        
        return [mean, median, stddev]

    def plot(self, color='r',lw=1,ls='solid'):
        
        for ap in self.ap:
            # Pintamos las aperturas
            ap.plot(color=color,ls=ls,lw=lw)
    
    def hist(self, bins=None, range=None, density=False, weights=None, 
             cumulative=False, bottom=None, histtype='bar', align='mid', 
             orientation='vertical', rwidth=None, log=False, color=None,
             label=None, stacked=False):
        
        import matplotlib.pyplot as plt
        
        plt.hist(self.counts, bins=bins, range=range, density=density, weights=weights, 
             cumulative=cumulative, bottom=bottom, histtype=histtype, align=align, 
             orientation=orientation, rwidth=rwidth, log=log, color=color,
             label=label, stacked=stacked)



# =============================================================================
# Algoritmos de centrado 
# =============================================================================

def Centroid(data,coords,r=3,method='1dg',mask=None):
    """
    Entrada: imagen, lista de coordenadas. Busca el centroide en una region de radio "r" entorno a la posicion dada en "coords".
    Salida: array con las posisiones ajustadas y distancia entre las posisiones "d"

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
        
    coords : numpy array
        n coordinates in a (n,2) numpy array.
    
    r : number, optional
        Size in pixels of the section arround the start where perform the centroid. The default is 3.
    
    method : photutils.centroid
        com: Calculates the object “center of mass” from 2D image moments.
        quadratic: Calculates the centroid by fitting a 2D quadratic polynomial to the data.
        1dg: Calculates the centroid by fitting 1D Gaussians to the marginal x and y distributions of the data.
        2dg: Calculates the centroid by fitting a 2D Gaussian to the 2D distribution of the data.

    Returns
    -------
    None.

    """    
    
    # =============================================================================
    # Paquetes utilizados     
    # =============================================================================
    from photutils.centroids import centroid_com, centroid_quadratic,centroid_1dg, centroid_2dg
    from astropy.stats import sigma_clipped_stats
    import numpy as np
    
    
    # Diccionario con los distintos metodos de centrados     
    cent = {'com': centroid_com, 'quadratic': centroid_quadratic, '1dg': centroid_1dg, '2dg':centroid_2dg}
    # Posiciones iniciales
    
    px = coords[0]
    py = coords[1]

    # Vamos a definir una seccion de los datos
    xmas = np.int(px+r)
    xmenos = np.int(px-r)
    ymas = np.int(py+r)
    ymenos = np.int(py-r)
    
    sec = data[ymenos:ymas,xmenos:xmas]

    #Calculamos el cielo solo dentro de la dregion selectionada y lo restamos        
    median = sigma_clipped_stats(sec, sigma=3.0)[1]
    sec = sec - median
    
    x_s, y_s = cent[method](sec, mask=mask)
    
    x = x_s + px - r
    y = y_s + py - r
    
    fit_coords = [x,y]
    
    return fit_coords


# =============================================================================
# Photometry along a line
# =============================================================================


class LinePhot:
    def __init__(self,x_a, x_b, y_a, y_b, xmin, xmax, w, h):

        # recta que pasa por los dos puntos dados
        import numpy as np
        m = (y_b-y_a)/(x_b-x_a)
        y0 = -(y_b-y_a)/(x_b-x_a)*x_a + y_a
        ang = np.arctan(m)
    
        # Aperturas
        from photutils import RectangularAperture

        delta_x =np.cos(ang)
        delta_y =np.sin(ang)
        
        # Loop para las aperturas a lo largo de la recta
        x = xmin
        y = xmin * m + y0
        
        self.ap =[]
        
        while x < xmax:
            self.ap.append(RectangularAperture([x, y], w, h, theta=ang))
            x = x + delta_x*w
            y = y + delta_y*w
         
        # Guardamos la variables como globales para usarlas en la ultima definicion
        self.m = m
        self.y0 = y0
        self.ang = ang

    def plot(self,color=None,ls='solid'):
        # =============================================================================
        # Paleta de colores
        # =============================================================================
        if color==None:
            from matplotlib.pyplot import cm
            paleta=cm.rainbow(np.linspace(0,1,len(self.ap)))
        else:
            paleta=[color]*(len(self.ap))
        
        for i in self.ap: 
            i.plot(color=paleta[self.ap.index(i)],ls=ls)

    
    def phot(self, data, method,subpixels=None,sigma=3):
        from astropy.stats import sigma_clipped_stats        
        
        counts = []
        x = []
        for i in self.ap: 
            sec_data = SecImData(i,data,method,subpixels=subpixels)
            # media sigma clip
            median = sigma_clipped_stats(sec_data, sigma=sigma)[1] 
            counts.append(median)
            x.append(i.positions[0])
        
        return np.array(x), np.array(counts)
    
    def centerLine(self,x):
        return self.m * x + self.y0
    
    def bondary(self,x,px,py):
        m = -1/self.m
        y0 = py - px*m
        return x*m +y0
        
    
    
# =============================================================================
# Sersic + exponencial
# =============================================================================


class SersExp:
    
    def Sersic(self,amplitude,r_eff,n):
        from astropy.modeling.models import Sersic1D
        self.Sers = Sersic1D(amplitude=amplitude, r_eff=r_eff, n=n)
        
    def SersFixed(self,amplitude=False,r_eff=False,n=False):
        self.Sers.amplitude.fixed = amplitude
        self.Sers.r_eff.fixed = r_eff
        self.Sers.n.fixed = n
        
    def SersBound(self,amplitude=(None,None),r_eff=(None,None),n=(None,None)):
        self.Sers.amplitude.bounds = amplitude
        self.Sers.r_eff.bounds = r_eff
        self.Sers.n.bounds = n

    
    def Exponential(self,amplitude,tau):
        from astropy.modeling.models import Exponential1D
        self.Exp = Exponential1D(amplitude=amplitude, tau=tau)
        
    def ExpFixed(self,amplitude=False, tau=False):
        self.Exp.amplitude.fixed = amplitude
        self.Exp.tau.fixed = tau
        
    def ExpBound(self,amplitude=(None,None),tau=(None,None)):
        self.Exp.amplitude.bounds = amplitude
        self.Exp.tau.bounds = tau

    
    # def Fitter(self,fitter = None):
    #     if fitter == None:
    #         from astropy.modeling.fitting import SLSQPLSQFitter 
    #         self.fitter = SLSQPLSQFitter()
    #     else:
    #         from astropy.modeling.fitting import fitter 
    #         self.fitter = fitter()

    def Fit(self,r,ct,weights = None, sigma = 3):
        from astropy.stats import sigma_clip
        from astropy.modeling.fitting import FittingWithOutlierRemoval
        from astropy.modeling.fitting import SLSQPLSQFitter 

        init = self.Sers + self.Exp  
        fitter = SLSQPLSQFitter()
        fit = FittingWithOutlierRemoval(fitter, sigma_clip, niter=3, sigma=sigma)

        return fit(init,r,ct,weights=weights)

    
# =============================================================================
# HI surface brightness (surf. mass density) profiles form VLA data
# =============================================================================

def surfDensVLAHI(data, bmaj, bmin, freq):
    '''
    Function to get the surface mass density of HI from VLA intensity maps in Jy/beam km/s
    
    First converts between surface brightness (flux density per area) and brightness temperature.
    Jy/beam to K by specifying the beam area
    https://docs.astropy.org/en/stable/units/equivalencies.html
    https://docs.astropy.org/en/stable/api/astropy.units.equivalencies.brightness_temperature.html#astropy.units.equivalencies.brightness_temperature
    
    Then converts k*km/s to Msun/pc^2 

    Parameters (from image header)
    ----------
    data : astropy.table.column.Column
        Column with the HI from VLA intensity maps in Jy/beam km/s
    bmaj : float
        BMAJ in deg
    bmin : float
        BMIN in deg
    freq : float
        RESTFREQ in Hz

    Returns
    -------
    Surface mass density in Msun/pc^2 as astropy.table.column.Column format

    '''
    
    import numpy as np
    from astropy import units as u
    
    bmaj = bmaj.to(u.arcsec)
    bmin = bmin.to(u.arcsec)
    fwhm_to_sigma = 1./(8*np.log(2))**0.5
    beam_area = 2.*np.pi*(bmaj*bmin*fwhm_to_sigma**2)
    '''
    freq = freq.to(u.GHz)
    equiv = u.brightness_temperature(freq)
    
    brigTemp = (data*u.Jy/beam_area).to(u.K, equivalencies=equiv)
    
    n_HI = 1.823e18 * brigTemp/u.K # HI atoms/cm^2  N(HI)[cm-2] =􏰂 1.82􏰉x10ˆ18 I_HI [K km/s ].
    
    atomHI = 1.673534e-27*u.kg # HI is neutral hidrogen = 1 p+ and 1 n-; same as regular H
    
    
    surfMassDens = (n_HI*atomHI/(u.cm*u.cm))
    surfMassDens = surfMassDens.to(u.Msun/(u.parsec*u.parsec))
    '''
    
    # Formula alternativa (eq. 2 en https://arxiv.org/pdf/2110.01618.pdf)
    # Da el mismo resultado de lo de arriba SIN corregir de inclinación 
    # (realmente yo tampoco corrijo arriba de inclinacion asi que tiene sentido)
    surfMassDens = 8794 * data /(bmaj*bmin)    
    
    #Si corrigiera de inclinacion: surfMassDens = 8794 * data * np.cos(88.5* np.pi / 180. )/(bmaj*bmin)


    return surfMassDens     

    
# =============================================================================
# Look for GAMA groups within a RA, DEC region
# =============================================================================

def table_gama_reg(ramin, ramax, decmin, decmax, tablepath):
    '''
    Gets GAMA groups within a given region in the sky.

    Parameters
    ----------
    ramin : float
        RA coord. min. value of the region in deg.
    ramax : float
        RA coord. max. value of the region in deg.
    decmin : float
        DEC coord. min. value of the region in deg.
    decmax : float
        DEC coord. max. value of the region in deg.
    tablepath : str
        Path of the GAMA table with the groups where to extract the targets.

    Returns
    -------
    None.

    '''



# =============================================================================
# Circular profiles
# =============================================================================

class CircularApProfile:
     """ Descripcion
     """
     def __init__(self, image, Xc, Yc, nbins, r_min, npix):        
        
         self.image=image
         self.Xc = Xc
         self.Yc = Yc
         self.nbins = nbins
         self.npix = npix
         self.r_min = r_min 
        
         #=======================================================================================================
         #  Empezamos definiendo los parametros para las aperturas
         #=======================================================================================================
         import numpy as np
         i=np.log10(r_min)
         l=np.log10(self.npix)
       
         step= np.logspace(i,l,self.nbins)
       
         self.r=step
         
         #==============================================================================
         # Aperturas
         #==============================================================================
       
         # La primera apertura es una ellipse y el resto son anillos concentricos
         #self.ap=[CircularAperture((self.Xc,self.Yc), r=self.r[0])]
         self.ap=[]
        
         for i in range(len(self.r)-1):
             self.ap.append(CircularAnnulus((self.Xc,self.Yc), r_in=self.r[i], r_out=self.r[i+1]))

         #=========================================================================================================

     def plot(self,image=None, color=None, ls='solid', cmap='viridis', vmin=None, vmax=None):
         # Paleta de colores
         if color==None:
             from matplotlib.pyplot import cm
             paleta=cm.rainbow(np.linspace(0,1,len(self.r)))
         else:
             paleta=[color]*(len(self.ap))
            
         for i in range(len(self.ap)):
             self.ap[i].plot(color=paleta[i],ls=ls)
         #=========================================================================================================xw
     
         
     def sbprofile(self, zp, pixscale, method='center',rms_sky=0, A=0, dimming=0, k_corr=0, subpixels=None, sigma=3):
         from astropy.stats import sigma_clipped_stats
         # =============================================================================
         # Primero calculamos la mediana dentro de cada anillo         
         # =============================================================================
         cts = []
         c_err = []
         for i in self.ap:
             # seccion dentro de la apertura
             sec_data = SecImData(i, self.image, method, subpixels=subpixels)
             # media sigma clip
             median = sigma_clipped_stats(sec_data, sigma=sigma)[1]  
             cts.append(median)

             # =============================================================================
             # Errores         
             # =============================================================================        
             # error Poissoniano
             error_poiss=np.sqrt(np.nanmean(np.square(median - sec_data))/np.size(sec_data))      
             # error bin
             error_bin = rms_sky/np.sqrt(float(np.size(sec_data)))
                 
             # Cuentas
             # error total
             err_t = np.sqrt(np.square(error_poiss)+np.square(error_bin))
             c_err.append(err_t)
         
         c_err = np.array(c_err) 
         cts = np.array(cts) 
 
         # Magnitudes
         # error total
         up = cts + c_err
         down = cts - c_err
         
         # Pasamos a magnitudes
         error_down = -2.5*np.log10(cts/up) 
         error_up = -2.5*np.log10(down/cts)



         # Definimos "r" como el punto medio de la anchura de la apertura
         c2sb = Counts2Sb(zp, pixscale, A=0, dimming=0, k_corr=0)
         r = self.r[:-1]+np.diff(self.r)/2
         sb = c2sb.do(cts)
         
         vec = [r, sb, error_up, error_down, cts, c_err]
         names = ['r', 'sb','err up','err down','ct','ct err']

        
         return Table(vec, names=names, meta={'Surface Brightnes': 'table info'})


        
     
     def model2d(self, method, subpixels=None, sigma=3):
         from astropy.stats import sigma_clipped_stats
         # =============================================================================
         # Primero calculamos la mediana dentro de cada anillo         
         # =============================================================================
         cts = []
         
         for i in self.ap:
             # seccion dentro de la apertura
             sec_data = SecImData(i, self.image, method, subpixels=subpixels)
             # media sigma clip
             median = sigma_clipped_stats(sec_data, sigma=sigma)[1]  
             cts.append(median)

         # Definimos "r" como el punto medio de la anchura de la apertura
            
         r = np.insert(self.r[1:]+np.diff(self.r)/2,0,self.r[0]/2)
             
         # =============================================================================
         # Luego costruimos el modelo 2D         
         # =============================================================================
         from scipy.interpolate import LSQUnivariateSpline   
         from photutils.isophote import EllipseGeometry
         # the target grid is spaced in 0.1 pixel intervals so as
         # to ensure no gaps will result on the output array.
         finely_spaced_sma = np.arange(r[0], r[-1], 0.1)

         # interpolate ellipse parameters
         # End points must be discarded, but how many?
         # This seems to work so far
         nodes = r[2:-2]
         shape = self.image.shape
         
         intens_array = LSQUnivariateSpline(r, cts, nodes)(finely_spaced_sma)
         x0 = self.Xc
         y0 = self.Yc

         result = np.zeros(shape=shape)
         weight = np.zeros(shape=shape)


         # for each interpolated ring, generate intensity values on the
         # output image array
    
         for index in range(1, len(finely_spaced_sma)):
             sma0 = finely_spaced_sma[index]
             geometry = EllipseGeometry(x0, y0, sma0, eps=0, pa=0)
    
             intens = intens_array[index]
    
             # scan angles. Need to go a bit beyond full circle to ensure
             # full coverage.
             r = sma0
             phi = 0.
             while phi <= 2*np.pi + 0.05:
                 # we might want to add the third and fourth harmonics
                 # to the basic isophotal intensity.
                 harm = 0.
    
                 # get image coordinates of (r, phi) pixel
                 x = r * np.cos(phi) + x0
                 y = r * np.sin(phi) + y0
                 i = int(x)
                 j = int(y)
    
                 if (i > 0 and i < shape[1] - 1 and j > 0 and j < shape[0] - 1):
                     # get fractional deviations relative to target array
                     fx = x - float(i)
                     fy = y - float(j)
    
                     # add up the isophote contribution to the overlapping pixels
                     result[j, i] += (intens + harm) * (1. - fy) * (1. - fx)
                     result[j, i + 1] += (intens + harm) * (1. - fy) * fx
                     result[j + 1, i] += (intens + harm) * fy * (1. - fx)
                     result[j + 1, i + 1] += (intens + harm) * fy * fx
    
                     # add up the fractional area contribution to the
                     # overlapping pixels
                     weight[j, i] += (1. - fy) * (1. - fx)
                     weight[j, i + 1] += (1. - fy) * fx
                     weight[j + 1, i] += fy * (1. - fx)
                     weight[j + 1, i + 1] += fy * fx
    
                     # step towards next pixel on ellipse
                     phi = max((phi + 0.75 / r), 0.05)
                     r = max(geometry.radius(phi), 0.5)
                 # if outside image boundaries, ignore.
                 else:
                     phi = max((phi + 0.75 / r), 0.05)
                     r = max(geometry.radius(phi), 0.5)

    
         # zero weight values must be set to 1.
         weight[np.where(weight <= 0.)] = 1.
    
         # normalize
         result /= weight
 
         return result



#********************************************************************************
# ELLIPTICAL APERTURES RADIAL PROFILE 
#********************************************************************************
class EllipticalApProfile:
    """ Descripcion
    """
    def __init__(self, image, Xc, Yc, sma_arr, ell_arr, pa_arr):
        
        
        self.image=image
        self.Xc = Xc
        self.Yc = Yc
        self.ell = ell_arr
        self.theta = pa_arr
        self.a = sma_arr
        
        #=======================================================================================================
        #  Empezamos definiendo los parametros para las aperturas
        #=======================================================================================================      
        self.b=(1-self.ell)*self.a

        #==============================================================================
        # Aperturas
        #==============================================================================
        from photutils import EllipticalAperture
        from photutils import EllipticalAnnulus
        
        # La primera apertura es una ellipse y el resto son anillos concentricos
        self.ap=[EllipticalAperture((self.Xc,self.Yc), a=self.a[0], b=self.b[0],theta=self.theta[0])]
        
        for i in range(len(self.a)-1):
            self.ap.append(EllipticalAnnulus((self.Xc,self.Yc), a_in=self.a[i],a_out=self.a[i+1], b_out=self.b[i+1],theta=self.theta[i]))

        #=========================================================================================================

    def plot(self,image=None,color=None,ls='solid',cmap='viridis',vmin=None,vmax=None):
        import matplotlib.pyplot as plt
        # Si no le damos la imagen, pinta la imagen original del analisis.
        if  np.shape(image)==() :
            image=self.image

        # Pintamos la imagen
        plt.imshow(image,cmap=cmap,vmin=vmin,vmax=vmax,origin='lower')
        

        # Paleta de colores
        if color==None:
            from matplotlib.pyplot import cm
            paleta=cm.rainbow(np.linspace(0,1,len(self.a)))
        else:
            paleta=[color]*(len(self.a))
            
        for i in range(len(self.ap)):
            self.ap[i].plot(color=paleta[i],ls=ls)
        
        #=========================================================================================================
            
    def surfbrigth(self, zp, pixscale, method='center',rms_sky=0, A=0, dimming=0, k_corr=0, subpixels=None, sigma=3):
        from astropy.table import Table
        from astropy.stats import sigma_clipped_stats
        
        self.zp = zp
        self.rms_sky = rms_sky
        self.A = A
        self.pixelscale = pixscale
        self.dimming = dimming
        self.k_corr = k_corr
        #==================================================================================
        #  El primer elemento es la elipse central. El resto de las aperturas son anillos
        #==================================================================================
        cts = []
        c_err = []
        sb = []
        error_up = []
        error_down = []
        noise = []

        for i in range(len(self.ap)):
            mask=self.ap[i].to_mask(method=method,subpixels=subpixels)
            mask_data=mask.data
            sec=mask.cutout(self.image)
            sec_weight=sec*mask_data
            sec_data=sec_weight[mask_data>0]            
            
            # =============================================================================
            # Intensity of the bin
            # =============================================================================
            median= sigma_clipped_stats(sec_data, sigma=sigma)[1]
            cts.append(median)

            # =============================================================================
            # Errores
            # =============================================================================
            # error Poissoniano
            error_poiss= np.sqrt(np.nanmean(np.square(median - sec_data))/float(np.size(sec_data)))      
            # error bin
            error_bin=self.rms_sky/np.sqrt(float(np.size(sec_data)))
            noise.append(error_bin)
        
            # error total
            err_t = np.sqrt(np.square(error_poiss)+np.square(error_bin))
            c_err.append(err_t)

            
        #==================================================================================
        
        c_err = np.array(c_err) 
        cts = np.array(cts) 

        # Magnitudes
        # error total
        up = cts + c_err
        down = cts - c_err
        
        # Pasamos a magnitudes
        error_down = -2.5*np.log10(cts/up) 
        error_up = -2.5*np.log10(down/cts)

        
        # Definimos "r" como el punto medio de la anchura de la apertura
        r=np.insert(self.a[1:]+np.diff(self.a)/2,0,self.a[0]/2)
        #c2sb = Counts2Sb(self.zp, self.pixelscale, A=0, dimming=0, k_corr=0)
        c2sb = Counts2Sb(self.zp, self.pixelscale, A=self.A, dimming=self.dimming, k_corr=self.k_corr)
        sb = c2sb.do(cts)
        
        # Senal a Ruido
        SNR = cts/noise
        
        
        vec = [r, sb, error_up, error_down, cts, c_err, SNR]
        names = ['r', 'sb','err up','err down','ct','ct err','SNR']
           
        return Table(vec, names=names, meta={'Surface Brightnes': 'table info'})
    
      
         
      
# =============================================================================
# Scale factor for subtracting a PSF scattered field of stars
# =============================================================================

class ScalePSF:
    """ Description
    """
    def __init__(self, image, psf, star_center, f_min,f_max, r_min = 50,r_max=400, psf_mask=None, psf_center ='mid', step=1):
        self.image=image
        self.psf=psf
        self.f_min = f_min
        self.f_max = f_max
        
        self.star_center = star_center
        
        # centro de la PSF
        if psf_center =='mid':
                psf_center= [int(psf.shape[0]/2)+1, int(psf.shape[1]/2)+1]
        
        #=======================================================================================================
        #  Empezamos definiendo los parametros para las aperturas
        #=======================================================================================================
        i=np.log10(r_min)
        l=np.log10(r_max+step)

        bins = int((r_max+step -r_min)/step)

        r = np.logspace(i,l,bins)

        # =============================================================================
        # Aperturas         
        # =============================================================================

        self.img_ap = []
        self.psf_ap = []

        for i in range(len(r)-1):
            self.img_ap.append(CircularAnnulus(star_center, r_in=r[i], r_out=r[i+1]))
            self.psf_ap.append(CircularAnnulus(psf_center, r_in=r[i], r_out=r[i+1]))

        # =============================================================================
        # Perfiles                     
        # =============================================================================
        img_profile = []
        psf_profile = []
        
        
        psf_masked = np.copy(psf)
        psf_masked[psf_mask!=0]=np.nan
            
        for i, j in zip(self.img_ap,self.psf_ap):
            # Image
            # seccion dentro de la apertura
            sec_data = SecImData(i, self.image, method='center')
            # media sigma clip
            median = sigma_clipped_stats(sec_data, sigma=3)[1]  
            img_profile.append(median)

            # PSF
            # seccion dentro de la apertura
            sec_data = SecImData(j, psf_masked, method='center')
            # media sigma clip
            median = sigma_clipped_stats(sec_data, sigma=3)[1]  
            psf_profile.append(median)
        
        self.r = r[:-1]+np.diff(r)/2
        self.img_profile = np.array(img_profile)
        self.psf_profile = np.array(psf_profile)
        
        self.selected = (self.img_profile < self.f_max) & (self.img_profile > self.f_min)

    
    def model_fluxRange(self,model=None):
        
        # =============================================================================
        # Factor de escala a partir del cociente de los perfiles en un rango dado
        # =============================================================================        
        
        self.scale_factor_dist = np.log10(self.img_profile[self.selected] / self.psf_profile[self.selected])
        self.scale_factor = sigma_clipped_stats(self.scale_factor_dist, cenfunc = 'median', sigma = 2.5)[0]  
    
        if model==None:
            model = np.zeros(self.image.shape)
            
        for i_psf in np.arange(self.psf.shape[0]):
            for j_psf in np.arange(self.psf.shape[1]):

                i = i_psf + int(self.star_center[1])- int(self.psf.shape[1]/2)
                j = j_psf + int(self.star_center[0])- int(self.psf.shape[0]/2)
                
                if (i > 0 and i < self.image.shape[1] - 1 and j > 0 and j < self.image.shape[0] - 1 and not np.isnan(self.psf[i_psf,j_psf])):
                    model[i,j]=self.psf[i_psf,j_psf] * 10**self.scale_factor

        return(model)

    
    
    def model_optimize(self,scale_factor=0, n_bs=100,model=None):

        # =============================================================================
        # Optimizacion del factor de escala
        # =============================================================================
        
        self.scale_factor_dist = []
        
        if scale_factor==0:
            for _ in range(n_bs):
                
                l = len(self.img_profile[self.selected])
                s = int(l*0.68)
                
                index = np.random.randint(0,l,size=s)
                
                def objective_function(a):
                    img_profile = self.img_profile[self.selected][index]
                    psf_profile = self.psf_profile[self.selected][index]
                    
                    return np.nansum(np.square((img_profile-a*psf_profile)*self.r[self.selected][index]))
             
                res = minimize_scalar(objective_function)
                self.scale_factor_dist.append(res.x)
                
            self.scale_factor = np.median(self.scale_factor_dist)
            print('Optimized sf :'+str(self.scale_factor))
        
        if scale_factor!=0:
            self.scale_factor = scale_factor
            print('Manual sf :'+str(self.scale_factor))


        if model==None:
            model = np.zeros(self.image.shape)
        # =============================================================================
        # Modelo 2D         
        # =============================================================================
        
        for i_psf in np.arange(self.psf.shape[0]):
            for j_psf in np.arange(self.psf.shape[1]):

                i = i_psf + int(self.star_center[1])- int(self.psf.shape[1]/2)
                j = j_psf + int(self.star_center[0])- int(self.psf.shape[0]/2)
                
                if (i > 0 and i < self.image.shape[1] - 1 and j > 0 and j < self.image.shape[0] - 1 and not np.isnan(self.psf[i_psf,j_psf])):
                    model[i,j]=self.psf[i_psf,j_psf] * self.scale_factor

        return(model)


    def plot(self,savefig=None):
        
        fig, ax = plt.subplots(2,2,figsize=(6,6))
        
        # =============================================================================
        # star         
        # =============================================================================
        coutout = self.img_ap[-1].to_mask(method='exact').cutout(self.image)
        # aperturas centradas en el corte de la imagen
        ax[0][0].imshow(coutout, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=10))
        x_c ,y_c= coutout.shape[0]/2,coutout.shape[1]/2

        # Pintamos el anillo interior y el exterior
        # Interior
        circ_int = CircularAperture([x_c,y_c], self.img_ap[0].r_in)
        circ_int.plot(axes=ax[0][0],ls='--',color='r')
    
        # Exterior
        circ_int = CircularAperture([x_c,y_c], self.img_ap[-1].r_out)
        circ_int.plot(axes=ax[0][0],ls='--',color='r')
        # Quitamos los ticks
        ax[0][0].set_xticks([])
        ax[0][0].set_yticks([])
        
        # =============================================================================
        # PSF         
        # =============================================================================
        coutout = self.psf_ap[-1].to_mask(method='exact').cutout(self.psf)
        # aperturas centradas en el corte de la imagen
        ax[0][1].imshow(coutout, origin='lower', cmap='viridis',norm=LogNorm(vmin=0,vmax=0.00001))
        x_c ,y_c= coutout.shape[0]/2,coutout.shape[1]/2
        # Pintamos el anillo interior y el exterior
        # Interior
        circ_int = CircularAperture([x_c,y_c], self.psf_ap[0].r_in)
        circ_int.plot(axes=ax[0][1],ls='--',color='r')
    
        # Exterior
        circ_int = CircularAperture([x_c,y_c], self.psf_ap[-1].r_out)
        circ_int.plot(axes=ax[0][1],ls='--',color='r')
        
        # Quitamos los ticks
        ax[0][1].set_xticks([])
        ax[0][1].set_yticks([])
        
        # =============================================================================
        # Perfil                 
        # =============================================================================
        ax[1][0].plot(self.r,self.img_profile,'-o', markersize=2, label='Star profile')
        ax[1][0].plot(self.r,self.psf_profile*self.scale_factor,'--o',markersize=2,label='PSF profile')

        #ax[1][0].scatter(self.r[self.selected],self.img_profile[self.selected],s=2)
        #ax[1][0].scatter(self.r[self.selected],self.psf_profile[self.selected]*10**self.scale_factor,s=2)


        ax[1][0].axhline(self.f_max,ls='--',color='k')
        ax[1][0].axhline(self.f_min,ls='--',color='k')
        ax[1][0].legend()
        ax[1][0].set_xscale('log')
        ax[1][0].set_yscale('log')
        ax[1][0].set_xlim([1, self.img_ap[-1].r_out+30])

       
        # =============================================================================
        # Histograma factor de escala                 
        # =============================================================================
        ax[1][1].hist(self.scale_factor_dist,bins='auto')
        ax[1][1].axvline(self.scale_factor,color='r',ls='--')
        
        
        plt.tight_layout()
        plt.show()  

        if savefig!=None:
            plt.savefig(savefig)





# =============================================================================
# Function: returns distances map from the centre of each galaxy [Mireia]
# =============================================================================

def dist_ellipse_map(N, xc, yc, ratio, pos_ang):  # MIREIA's!!
    import numpy as np
    #Convert POS_ANG to radians
    ang = pos_ang * np.pi / 180.
    cosang = np.cos(ang)
    sinang = np.sin(ang)
    if len(N) == 2 :
        nx = N[0]
        ny = N[1]
    elif len(N)==1 :
        ny = N[0]
        nx = N[0]
    x = np.arange(0, nx,1) - [xc-1.]*nx
    y = np.arange(0, ny,1) - [yc-1.]*ny
    im = np.zeros([ny, nx])
    # Rotate pixels to match ellipse orientation
    xcosang = x*cosang
    xsinang = x*sinang
    for ii in np.arange(0, ny, 1):
         xtemp =  xcosang + y[ii]*sinang
         ytemp = -xsinang + y[ii]*cosang
         # ratio = 1 - elll; ell = 1-b/a
         # eq. ellipse: 1 = (x/a)**2 + (y/b)**2
         im[ii,:] = np.sqrt((ytemp/ratio)**2 + xtemp**2 ) # SMA
    return im



def dist_galaxies(image, x_gals, y_gals, ratio_gals, pa_gals):

    ny,nx = np.shape(image)
    dist = []
    for tt in np.arange(0, len(x_gals)):
        grid = dist_ellipse_map([nx,ny], x_gals[tt], y_gals[tt], ratio_gals[tt], pa_gals[tt])
        grid = np.array(grid, dtype = int)
        dist.append(grid)

    r = np.amin(dist, axis = 0)
    return r #* pixscale * arcsec_to_kpc[ii] # 2D array of min dist r to the com of gals 



# =============================================================================
# Function: extracts the (counts) profile using the distances map
# =============================================================================

def dist_ellipse_prof(masked_data, dist_map, sma0, step, nbins, SkyRMS, zp, pixscale, A=0, dimming=0, k_corr=0):
    '''
    

    Parameters
    ----------
    masked_data : np array 
        Group image values from the masked_data fits file.
    dist_map : np array 
        dist_map image values from the dist_map fits file.
    sma0 : int
        The starting value for the first apperture length (pixels).
    step : float
        The step value for growing/shrinking the apperture (pixels). 
    nbins : int
        Number of appertures for the profile.
    SkyRMS : float
        Sky background value.
    zp : float
        Zero point [mag]
    pixscale : float
        Pixel scale in arcsec/pix.
    A : float, optional
        Galactic absorption or extinction for a given band. The default is 0.
    dimming : float, optional
        Dimming: 10*np.log10(1+z) [mag/arcsec^2]. 
        Because the diming in linear scale is (1+z)**-4 [Calvi+2014, doi: 10.1088/0004-637X/796/2/102]. 
        The default is 0.
    k_corr : float, optional
        K correction [J. Loveday et al. (2012)]. When obs. redshifted galaxies there is an spectral shift, this is a change of the wavelength values as lambda*(1+z) so as wavelenth range. 
        The default is 0.

    Returns
    -------
    astropy.table
        Table with the surface brightness values and errors.

    '''        
    # Function to get a table with the profile values + errors    

    from astropy.table import Table
    from astropy.stats import sigma_clipped_stats

    ell_bin = [sma0*(1. + step)]
    ct = []
    r = []
    c_err = []
    noise = []
    
    for nbin in range(nbins):
        # =============================================================================
        # Bins & intensity (ct)              
        # =============================================================================
        ell_bin.append(ell_bin[nbin]*(1. + step))
        data_bin = masked_data[(dist_map >= ell_bin[nbin]) & (dist_map < ell_bin[nbin+1])]
        
        median = sigma_clipped_stats(data_bin, sigma=1, mask_value=np.nan)[1]
        ct.append(median)
        
        r.append(ell_bin[nbin] + (ell_bin[nbin+1]- ell_bin[nbin])/2)
    

        # =============================================================================
        # Errors
        # =============================================================================
        # Poissonian error 
        error_poiss = np.sqrt(np.nanmean(np.square(median - data_bin))/float(np.size(data_bin)))      
        # Bin error
        error_bin = SkyRMS/np.sqrt(float(np.size(data_bin)))
        noise.append(error_bin)
    
        # Total error
        err_t = np.sqrt(np.square(error_poiss)+np.square(error_bin))
        c_err.append(err_t)

    #==================================================================================
    
    # Lists to arrays
    r = np.array(r)
    ct = np.array(ct) 
    ct_err = np.array(c_err) 

    # Signal to noise ratio
    SNR = ct/noise
 
    #==================================================================================
    # Convert counts to sb
    c2sb = Counts2Sb(zp, pixscale, A=A, dimming=dimming, k_corr=k_corr)
    sb = c2sb.do(ct)
    sb_err = c2sb.do_errors(ct_err)
    
    # Save Astropy.table        
    vec = [r, ct, ct_err, SNR, sb, sb_err[0], sb_err[1]]
    names = ['r', 'ct','ct err','SNR', 'sb','err up','err down']
       
    return Table(vec, names=names, meta={'Counts and SB Profile': 'table info'})
 
    
 
# =============================================================================
# Function: extracts radius at which the ICL dominates the flux from the distances map
# =============================================================================

def  findRadius(image, center = None):
     ## For the center = x, y
     imsize = np.shape(image)

     if center == None:
        center = [int(imsize[1]/2.), int(imsize[0]/2.)]

     distIm    = dist_ellipse_map(imsize,  center[0], center[1], 1, 0)
     d_vec     = np.arange(0, 2000, 10)
     profiles  = []
     radius    = []

     for nn in np.arange(0, len(d_vec)):
         ind = (distIm < d_vec[nn+1]) & (distIm > d_vec[nn])
         median_value = np.nanmedian(image[ind])

         if (len(image[ind][np.isnan(image[ind])]) > (len(image[ind])/2)):
             median_value = np.nan
         profiles.append(median_value)
         radius.append(d_vec[nn] + (d_vec[nn+1] - d_vec[nn])/2.)
         if ~np.isnan(median_value):
            break

     return radius[-1]
 

# =============================================================================
# Function: extracts the (counts) profile using the distances map
# =============================================================================

def isophotes_prof(masked_data, isophotes_im, Xc, Yc, SkyRMS, zp, pixscale, A=0, dimming=0, k_corr=0):
       
    # Function to get a table with the profile values + errors    

    from astropy.table import Table
    from astropy.stats import sigma_clipped_stats
    
    isoph_values = np.unique(isophotes_im)
    
    isoph_pixs_val = []
    isoph_pixs_coords = []
    ct = []
    r = []
    c_err = []
    noise = []
    
    
    # =============================================================================
    # For each ispophote:            
    # =============================================================================
    
    for i, isoph_value in enumerate(isoph_values):

        # 1. Calculate median counts value
        isoph_pixs_val.append(masked_data[isophotes_im==isoph_value])
        ct_median = sigma_clipped_stats(isoph_pixs_val[i], sigma=1, mask_value=np.nan)[1]
        ct.append(ct_median)
         
        ''' 
        # Opt 1. Calculate the median distance r [def as the median dist to all pix in a given isophote]
        isoph_pixs_coords.append(np.argwhere(isophotes_im==isoph_value))
        dist_pixs = []
         
        for c in isoph_pixs_coords[i]:
             d = np.sqrt( np.square(Xc-c[0]) + np.square(Yc-c[1]) )
             dist_pixs.append(d)                         
        r_median = sigma_clipped_stats(dist_pixs, sigma=3, mask_value=np.nan)[1]                            
        '''
        # Opt 2. The r distance is actually the isophote number
        r.append(isoph_value) 
    

        # =============================================================================
        # Errors
        # =============================================================================
        # Poissonian error 
        error_poiss = np.sqrt(np.nanmean(np.square(ct_median - isoph_pixs_val[i]))/float(np.size(isoph_pixs_val[i])))      
        # Bin error
        error_bin = SkyRMS/np.sqrt(float(np.size(isoph_pixs_val[i])))
        noise.append(error_bin)
    
        # Total error
        err_t = np.sqrt(np.square(error_poiss)+np.square(error_bin))
        c_err.append(err_t)

    #==================================================================================
    
    # Lists to arrays
    r = np.array(r)
    ct = np.array(ct) 
    ct_err = np.array(c_err) 

    # Signal to noise ratio
    SNR = ct/noise
 
    #==================================================================================
    # Convert counts to sb
    c2sb = Counts2Sb(zp, pixscale, A=A, dimming=dimming, k_corr=k_corr)
    sb = c2sb.do(ct)
    sb_err = c2sb.do_errors(ct_err)
    
    # Save Astropy.table        
    vec = [r, ct, ct_err, SNR, sb, sb_err[0], sb_err[1]]
    names = ['r', 'ct','ct err','SNR', 'sb','err up','err down']
       
    return Table(vec, names=names, meta={'Counts and SB Profile': 'table info'})
 
    
 
# =============================================================================
# Function to save images in .fits format   
# =============================================================================
 
def save_fits_image(data, header, save_path, new_file_name):
    """
    Function to write and save a .fits image

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    header : TYPE
        DESCRIPTION.
    save_path : TYPE
        DESCRIPTION.
    new_file_name : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    # Save distances map .fits
    hdu = fits.PrimaryHDU(data)  # we create a PrimaryHDU object to encapsulate the data:
    wcs = WCS(header)
    hdu.header.update(wcs.to_header())
    hdu.writeto(save_path+new_file_name, overwrite=True) # write to a new file




def flux_in_ellipse(data, Xc, Yc, maxsma, ell, pa):
    '''
    Estimate the amount of light (Flux in counts) within an elliptical aperture 

    Parameters
    ----------
    data : np.array
        Data where to measure the flux.
    Xc : float
        X coordinate [pix] of the centre of the ellipse as for Python.
    Yc : float
        X coordinate [pix] of the centre of the ellipse as for Python.
    maxsma : float
        Semi-major axis length of the ellipse [pix] .
    ell : float
        Ellipticity of the ellipse.
    pa : float
        Position angle of the ellipse.
    rms_sky : float, optional
        RMS of the sky in the image.

    Returns
    -------
    float
        Total flux within the ellipse in counts, associated error, and the 
        elliptical aperture

    '''
    b = (1-ell)*maxsma
    
    Ell_Ap = EllipticalAperture((Xc, Yc), a=maxsma, b=b, theta=ell)
    Ell_Section = SecImData(Ell_Ap, data, method='center', subpixels=None)
        
    # Sum all the values within the aperture to get the total flux
    flux = np.nansum(Ell_Section)
    # Error within the region due to the sky rms              
    #error = rms_sky/np.sqrt(float(np.size(Ell_Section)))
    error = np.sqrt(np.std(Ell_Section))

    return flux, error, Ell_Ap



def patch_mask_cold(im_chunk):
    '''
    Function to inject in the multiprocessing algorithm (parallelized 
    computing). Gnerates a chunk element with the detections in one patch
    
    NOT finished!! Needs tweaks!!

    Parameters
    ----------
    im_chunk : TYPE
        DESCRIPTION.

    Returns
    -------
    mask_chunk : TYPE
        DESCRIPTION.

    '''

    #threshold = detect_threshold(im_chunk, nsigma=1.1, sigclip_sigma=2, mask_value = np.nan)
    #mask_chunk = np.zeros_like(im_chunk)
    #try:
        # Segmentation - detecting sources
        #segm = detect_sources(im_chunk, threshold, npixels=1200,filter_kernel=kernel)
        # Deblending sources 
        #mask_chunk = deblend_sources(im_chunk, segm, npixels=1200, filter_kernel=kernel).data
    #except:
        #print('No object detected in patch ', patch)
    #return mask_chunk




def magnitudes(ct_array, zp, DMpc=0, z=0, k_corr=0, band='i'):
    '''
    Funtion to obtain the absolute magnitude of a galaxy/structure from a value in counts
    
        1. Get the aparent magnitude corrected by k_corrs
            m = -2.5*log10(ct_arr) + zp - k_corr
        2. Get the absolute magnitude using the distance to the object (in pc)
            M = m - 5*log10()

        abs_mag_sun= {'g':5.11, 'r':4.65, 'i':4.53} 
        #AB system from Table 3 in Willmer2018 https://iopscience.iop.org/article/10.3847/1538-4365/aabfdf


    Parameters
    ----------
    ct_array : np.array  
        Flux values in counts 
    zp : float
        Zero point of the image in mag
    scale : float
        Pixel sacle of the image in arcsec/pix
    DMpc : float (optional)
        Distance in Mpc to the object
    z : float (optional)
        Source redshift
    k_corr : float (optional)
        K-correction
    band : string (optional)
        i-band by default)
    Returns
    -------
    Dictionary of {apparent, absolute} magnitudes.

    '''            
    
    mag_array = -2.5*np.log10(ct_array) + zp - k_corr

    
    # If the source is nearby:
    if DMpc!=0:    
        # 
        abs_mag = mag_array-5*np.log10((DMpc*10**6)/10)
        
    
    # If source at some redshift, we should use the luminosity distabce method:
    if z!=0:
        from astropy.cosmology import FlatLambdaCDM
        import astropy.units as u
        
        cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)
        dist_lum = cosmo.luminosity_distance(0.2)
        
        abs_mag = mag_array-5*np.log10((dist_lum.value*10**6)/10)
        
    m = 'm_'+band
    M = 'M_'+band
    
    return {m: mag_array, M: abs_mag}



def luminosity(abs_mag, band='i'):
    '''
    Funtion to obtain the luminosity from the absolute magitude of a 
    galaxy/structure in a given band using the solar abs. mag as reference
           
        M = Msun - 2.5*log10(L/Lsun)
    
    ** Still have to figure out the luminosity errors
    
        Lum in [L_sun] -- this is the output units!
        L_SUN =  3.83e26 watts
        abs_mag_sun= {'g':5.11, 'r':4.65, 'i':4.53} 
        #AB system from Table 3 in Willmer2018 https://iopscience.iop.org/article/10.3847/1538-4365/aabfdf


    Parameters
    ----------
    ct_array : np.array  
        Flux values in counts
    band : string (optional)
        i-band by default)
    Returns
    -------
    Array of luminosity values in [L_sun].

    '''            
    abs_mag_sun = {'g':5.11, 'r':4.65, 'i':4.53} 

    luminosity = 10**(0.4*(abs_mag_sun[band]-abs_mag)) # this is in Lsun units
        
    #lum_err = ct_err_array*10**(0.4*(abs_mag_sun)) #in Lsun units No esta bien calculado
        
    return luminosity  # , lum_err




def massBell(luminosity, color, color_bands, band):
    '''
    
    Mass estimation for galaxies/groups/clusters/structures using the M/L and colors
    From Table 7 in Bell2003 https://iopscience.iop.org/article/10.1086/378847/pdf
    
    log10(M/L)= a_lambda + (b_lambda * color) ; M/L ratio is in solar units
    So then: 
        M = L * 10**(a_lambda + (b_lambda * color))
    
    
    Milky Way mass : 6.43 ± 0.63 x 10^10 Msun  (McMillan 2011)

    Parameters
    ----------
    luminosity : float
        Luminosity of the galaxy/structure in the corresponding band.
    color : float
        Color value (color = band1 - band2).
    color_bands : string
         Bands from which we have obtained the color i.e.: color_bands='g-r'  

    band : string
         Band in which we have obtained the luminosity i.e.: band='g'  

    Returns
    -------
    Mass value in M_sun units

    '''
    if color_bands=='g-r':
        a = {'g':-0.499, 'r':-0.306, 'i':-0.222}
        b = {'g':1.519, 'r':1.097, 'i':0.864}
        
    if color_bands=='g-i':
        a = {'g':-0.379, 'r':-0.220, 'i':-0.152}
        b = {'g':0.914, 'r':0.661, 'i':0.518}

    if color_bands=='r-i':
        a = {'g':-0.106, 'r':-0.022, 'i':0.006}
        b = {'g':1.982, 'r':1.431, 'i':1.114}

    return luminosity * 10**(a[band] + (b[band] * color))






def massTaylor(color_g_i, abs_mag_i):
    '''
    
    Mass estimation for galaxies/groups/clusters/structures using the (g-i) color 
    and the absolute magnitude in i-band
    From eq. 8 in Taylor+2011 (GAMA Stellar Masses paper):
        
    log10(Mass) = 1.15 + 0.70*(g-i)-0.4*M_i  
    
    Mass is in solar units
    M_i is the absolute magnitude in i-band
    
    Milky Way mass : 6.43 ± 0.63 x 10^10 Msun  (McMillan 2011)

    Parameters
    ----------
    color_g_i : float
        Color value (color = g - i).
    abs_mag_i : float
         Absolute magnitude of the galaxy/structure in i-band 

    Returns
    -------
    Mass value in M_sun units

    '''

    return 10**(1.15 + 0.70*color_g_i - 0.4*abs_mag_i)



def stellmass2halo(stellar_mass):
    '''
    Convertion from stellar mass to halo mass of a system following the recipe in Huang+20:
        log10(Mstell) [Msun] = 0.589 * (log10(Mhalo) - 13.5) + 11.844

    Parameters
    ----------
    stellar_mass : float
        Stellar mass of the group/cluster in M_sun units.

    Returns
    -------
    Halo (virial) mass of the system in M_sun units.

    '''
    halomass = 10**(((np.log10(stellar_mass) - 11.844)/0.589) + 13.5)
    
    return halomass




# =============================================================================
# Funcion para los dobles ejes
# =============================================================================

def kpc_arcsec_axis(ax,kpc_ticks, arcsec_ticks, arcsec_lim, size_label, size_ticks,scale, DMpc,arcsec_label=True,kpc_label=True, z=None, spatial_label= 'R'):
    
    from matplotlib.ticker import AutoMinorLocator


    kpc_ticks = np.array(kpc_ticks)
    arcsec_ticks = np.array(arcsec_ticks)
    # Kpc
    kpc_pos = Kpc2Arcsec(kpc_ticks, DMpc)
    
    # Imponemos que se pinten solo los ticks menores que el limite en arcsec
    valid_range = (kpc_pos > arcsec_lim[0]) & (kpc_pos < arcsec_lim[1])
    kpc_pos = kpc_pos[valid_range]
    kpc_ticks = kpc_ticks[valid_range]
    ax.set_xticks(kpc_pos)
    ax.set_xticklabels([])
    if kpc_label:
        ax.set_xlabel(spatial_label+" [kpc]",size = size_label)
        ax.set_xticklabels(kpc_ticks)

    ax.set_xlim(arcsec_lim)
    
    # arcsec    
    # Imponemos que se pinten solo los ticks menores que el limite en arcsec
    valid_range = (arcsec_ticks > arcsec_lim[0]) & (arcsec_ticks < arcsec_lim[1])
    arcsec_ticks = arcsec_ticks[valid_range]

    ax_t = ax.twiny()
    ax_t.tick_params(axis="both",direction="out")
    ax_t.tick_params(axis='x')
    ax_t.set_xlim(arcsec_lim)
    ax_t.set_xticks(arcsec_ticks)
    ax_t.set_xticklabels([])
    ax_t.xaxis.set_minor_locator(AutoMinorLocator())

    if arcsec_label:
        ax_t.set_xlabel(spatial_label+" [arcsec]", size = size_label)
        ax_t.set_xticklabels(arcsec_ticks)














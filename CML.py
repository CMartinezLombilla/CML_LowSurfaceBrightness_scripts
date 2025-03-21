#*****************************************************************************
#                             SB PROFILES FUNCTIONS 
#-----------------------------------------------------------------------------
#*****************************************************************************

# Packages versions
import numpy as np
import matplotlib.pyplot as plt 



# Versions
# numpy== 1.19.2
# matplotlib==3.3.2
# astropy==4.1
# photutils==1.0.1
# tkinter == 8.6

# conda install astropy==4.1
# conda install photutils==1.0.1


#==============================================================================
# Ignored warnings
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
            Zero point [mag/arcsec^2]
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
        
    def do(self,x):
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
        import numpy as np
        return - 2.5*np.log10(x) + self.zp - self.A - self.dimming - self.k_corr + 5*np.log10(self.pixscale)


    
def SecImData(ap,image,method,subpixels=None):
    '''
    Extracts a section (using any astropy.apperture shape) of a 2D array

    Parameters
    ----------
    ap : astropy.apperture class
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

    '''
    mask=ap.to_mask(method=method,subpixels=subpixels)
    mask_data=mask.data
    sec=mask.cutout(image)
    sec_weight=sec*mask_data
    sec_data = sec_weight[mask_data>0]
    # Quitamos los valores NaN
    sec_data = sec_data[np.isfinite(sec_data)]
    
    return sec_data



#==============================================================================
# Ajuste para buscar el centro        
#==============================================================================
class Center:

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
        
        from photutils import centroid_com #, centroid_1dg, centroid_2dg centroid_com
        
        #Las posiciones de los centroides ajustados
        x_s, y_s =centroid_com(self.sec)
        
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
        '''
        Main class definition to extract a surface brightness profile using rectangular appertures over an image. Here the appertures are built.     
        
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
        Class object with the appertures already defined.

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
        self.delta=delta
        self.orientation=orientation
        self.pixscale = pixscale
        
        #=======================================================================================================
        #  Empezamos definiendo los parametros para las aperturas
        #=======================================================================================================
        import numpy as np
        l=np.log10(self.npix)
        step=np.logspace(0,l,self.nbins)

        # Diccionario para todos los parametros sgun se avertical u horizontal.
        self.param={'horizontal':{'w' : np.diff(step), 'h' : np.ones(len(np.diff(step)))*self.height, 'x' : step[1:]-(np.diff(step))/2-step[0]-np.diff(step)[0]/2, 'y' : np.zeros(len(np.diff(step)))+self.delta, 'signo' : (-1,1)},            
               'vertical':{'w' : np.ones(len(np.diff(step)))*self.height,'h' : np.diff(step), 'x' : np.zeros(len(np.diff(step)))+self.delta, 'y' : step[1:]-(np.diff(step))/2-step[0]-np.diff(step)[0]/2, 'signo' : (1,-1), 'spatial' : ('z',step[1:]-(np.diff(step))/2-step[0]-np.diff(step)[0]/2)}}

        w=self.param[orientation]['w']
        h=self.param[orientation]['h']
        x=self.param[orientation]['x']
        y=self.param[orientation]['y']
        self.signo=self.param[orientation]['signo']

        #==============================================================================
        # Aperturas
        #==============================================================================
        from photutils import RectangularAperture

        self.ap_c=RectangularAperture([x[0]+self.Xc,y[0]+self.Yc], w=w[0], h=h[0],theta=0) # Apertura central
        self.ap_r=[] # Aperturas hacia la derecha
        self.ap_l=[] # Aperuras hacia la izquierda
        


        for i in range(len(x)-1):
            self.ap_r.append(RectangularAperture((x[i+1]+self.Xc,y[i+1]+self.Yc), w=w[i+1], h=h[i+1],theta=0))
            self.ap_l.append(RectangularAperture((self.signo[0]*x[i+1]+self.Xc,self.signo[1]*y[i+1]+self.Yc), w=w[i+1], h=h[i+1],theta=0))
        
        # Almancenamos las posiciones como variables globales y al diccionario para las tablas
        self.x=x+self.Xc
        self.y=y+self.Yc
        
        self.param['horizontal']={'spatial' : ('r', self.x)}
        self.param['vertical']={'spatial' : ('z', self.y)}


        #=========================================================================================================
        
    def plot(self,color=None,ls='solid',cent_ap=False):
        '''
        Plot of the profile appertures.

        Parameters
        ----------
        color : str, optional
            Color of the appeertures in the plot. The default is None.
        ls : str, optional
            Linestyle "solid", "dotted", "dashed", "dashdot", ... The default is 'solid'.
        cent_ap : bool, optional
            Plot a mark in the central position of the apperture. The default is False.

        Returns
        -------
        Plot of the appertures

        '''
        import matplotlib.pyplot as plt
                
        # Paleta de colores
        if color==None:
            from matplotlib.pyplot import cm
            paleta=cm.rainbow(np.linspace(0,1,len(self.x)))
        else:
            paleta=[color]*(len(self.x))
        
        # Pintamos las aperturas
        self.ap_c.plot(color=paleta[0],ls=ls)
        if cent_ap == True:   
            plt.text(self.x[0],self.y[0],0,ha='center',va='bottom')
            plt.scatter(self.x[0],self.y[0],marker='x',color='k',s=6)
                    
        for i in range(len(self.x)-1):
            self.ap_r[i].plot(color=paleta[i+1],ls=ls)
            self.ap_l[i].plot(color=paleta[i+1],ls=ls)
            
            if cent_ap == True:
                x_r=self.x[i+1]
                y_r=self.y[i+1]
                                
                # Numero
                plt.text(x_r,y_r,i+1,ha='center',va='bottom')
                # Una "x" en el centro
                plt.scatter(x_r,y_r,marker='x',color='k',s=6)
 
                
    def saveDS9(self, fname):
        '''
        Writes the appertures in a file as DS9 regions

        Parameters
        ----------
        fname : str
            File name.

        Returns
        -------
        A .reg file in the current folder.

        '''
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
                    
        for i in range(len(self.x)-1):
            # Aperturas derechas
            x,y = self.ap_r[i].positions
            w = self.ap_r[i].w
            h = self.ap_r[i].h
        
            l_r='box('+str(x)+','+str(y)+','+str(w)+','+str(h)+') # text={'+str(i+1)+'} dash=1'
            file.write(l_r+ "\n")
            
            # Aperturas izquierdas
            x,y = self.ap_l[i].positions
            w = self.ap_l[i].w
            h = self.ap_l[i].h
        
            l_l='box('+str(x)+','+str(y)+','+str(w)+','+str(h)+') # text={'+str(i+1)+'} dash=1'
            file.write(l_l+ "\n")
    
        file.close()
    
    
    
    def surfbrigth(self, method,subpixels=None,sigma=3):
        '''
        Extracts the surface brightness profile from the appertures defined in CentralProfile class.

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
        
        from astropy.table import Table
        from astropy.stats import sigma_clipped_stats        
        # Clase audiliar que define la transformacion de cuantas a brillo superficial         
        c2sb = Counts2Sb(self.zp,self.pixscale,self.A)

        #==================================================================================
        #  El primer elemento tanto a la derecha como a la izquierda es la apertura central
        #==================================================================================
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
        for i in range(len(self.ap_r)):
            # Aperturas derechas
            sec_data_r = SecImData(self.ap_r[i],self.image,method,subpixels=subpixels) 
            # media sigma clip
            median_r= sigma_clipped_stats(sec_data_r, sigma=sigma)[1]

            # Aperturas izquierdas
            sec_data_l = SecImData(self.ap_l[i],self.image,method,subpixels=subpixels)
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
        # Definimos ras para restar a la variable espacial y tenerla referida a 0
        ras = np.nanmin(self.param[self.orientation]['spatial'][1])
        spatial = self.param[self.orientation]['spatial'][1]
        spatial_name = self.param[self.orientation]['spatial'][0]
        
        vec = [spatial-ras, sb_r, sb_l,sb,error_up,error_down,c_r, c_l,c,c_err]
        names = [spatial_name, 'sb right', 'sb left','sb','err up','err down', 'ct right', 'ct left','ct','ct err']
        
        return Table(vec, names=names, meta={'Surface Brightnes': 'table info'})

#==============================================================================
# Perfiles radiales desplazados        
#==============================================================================
            
class ShiftedProfile:
    
    def __init__(self,image,Xc,Yc,nbins,npix,height,zp,pixscale,delta,rms_sky=0,A=0,orientation='horizontal'):
        
        '''
        Class definition to extract a surface brightness profile with an offset from the main axis of the galaxies using rectangular appertures over an image. 
        Here the appertures are built.    
        
        
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

    def plot(self,color=None,ls='solid',cent_ap=False):
        '''
        Plot of the profile appertures.

        Parameters
        ----------
        color : str, optional
            Color of the appeertures in the plot. The default is None.
        ls : str, optional
            Linestyle "solid", "dotted", "dashed", "dashdot", ... The default is 'solid'.
        cent_ap : bool, optional
            Plot a mark in the central position of the apperture. The default is False.

        Returns
        -------
        Plot of the appertures
        '''
        
        self.up.plot(color=color,ls=ls,cent_ap=cent_ap)
        self.down.plot(color=color,ls=ls,cent_ap=cent_ap)
    
    def saveDS9(self, fname):
        '''
        Writes the appertures in a file as DS9 regions

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
        Extracts the surface brightness profile from the appertures defined in ShiftedProfile class.

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

    def __init__(self, data, Xc, Yc, sma, eps, pa, astep=0.1, linear_growth=False, fix_center=False, fix_pa=False, fix_eps=False,ell_threshold=0.1):
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
        self.ellipse = Ellipse(data, self.geometry, threshold=ell_threshold)
        
    def find_center(self):
        '''
        Find the center coordinates of the galaxy.

        Returns
        -------
        If the algorithm is successful the (x, y) coordinates in this ~photutils.isophote.EllipseGeometry (i.e., the x0 and y0 attributes) instance will be modified.

        '''
        self.geometry.find_center(self.data)

    
    def fit_image(self, sma0=None, minsma=0.0, maxsma=None, step=0.1, conver=0.05, minit=10, maxit=50, fflag=0.7, maxgerr=0.5, sclip=3.0, nclip=0, integrmode='median', linear=None, maxrit=None, fix_center=False, fix_pa=False, fix_eps=False):
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

        # then we build a CentralEllipseSample instance and call the 
        # CentralEllipseFitter on it.
        sample = CentralEllipseSample(self.data, 0., geometry=g)
        fitter = CentralEllipseFitter(sample)
        center = fitter.fit()
        # Replace the new fit in the central region
        self.isolist[0] = center
        self.isolist.sort()

    def outer_regions(self,sma_cut,step):
        from photutils.isophote import EllipseGeometry
        from photutils.isophote import Ellipse

        # Indice de la isofota a partir de la cual ajustaremos denuevo
        ind = np.where(self.isolist.sma>=sma_cut)[0].min()

        # Los parametros iniciales seran los de la ultima elipse considerada
        x0 = self.isolist.x0[ind]
        y0 = self.isolist.y0[ind]
        sma = self.isolist.sma[ind]
        eps = self.isolist.eps[ind]
        pa = self.isolist.pa[ind]

        g = EllipseGeometry(x0, y0, sma, eps, pa)
       
        ellipse = Ellipse(self.data, geometry=g)

        outer_isolist = ellipse.fit_image(integrmode='median', step=step, minsma=sma,maxsma=self.maxsma, fflag=0.3, sclip=3.0, nclip=3, fix_center=True)

        # Agragamos el ajuste de la parte externa
        self.isolist = self.isolist[:ind] + outer_isolist
        
    def run_over(self,new_data, zp, pixscale, rms_sky=0, A=0, dimming=0, k_corr=0):
        from photutils.isophote import EllipseSample, Isophote, IsophoteList
        # from photutils.isophote.sample import CentralEllipseSample
        # from photutils.isophote.fitter import CentralEllipseFitter
        # from photutils.isophote import EllipseGeometry

        # =============================================================================
        # Bucle sonbre todas las isofotas 
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
        # Si no se especifica una isolist coge la del fir original
        if isolist == None:
            isolist = self.isolist
        # Clase audiliar que define la transformacion de cuantas a brillo superficial         
        c2sb = Counts2Sb( zp, pixscale, A, dimming, k_corr)
        
        # Leesmos las cuentas del ajuste
        sma = isolist.sma
        counts = isolist.intens
        counts_err = np.sqrt(np.square(isolist.rms))# + np.square(rms_sky))
        
        # Senal a Ruido
        SNR = isolist.intens/isolist.pix_stddev
        
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


        # Errores en magnitudes
        up = counts + counts_err
        down = counts - counts_err

        error_down = -2.5*(np.log10(counts/up))
        error_up = -2.5*(np.log10(down/counts))
        

        
        # =============================================================================
        # Reemplazamos los nan por 0 
        # =============================================================================
        error_up = np.nan_to_num(error_up)
        error_down = np.nan_to_num(error_down)

        
        # =============================================================================
        # Guardamos la tabla         
        # =============================================================================

        from astropy.table import Table
    
        names=['sb','err up','err down','ct','ct err','SNR','sma','xc', 'xc err', 'yc', 'yc err', 'ell', 'ell err', 'pa', 'pa err', 'b4', 'b4 err']
        vec = [sb, error_up, error_down, counts,counts_err,SNR,sma,xc, xc_err, yc, yc_err, ell, ell_err, pa, pa_err, b4, b4_err] 
    
        return Table(vec,names=names, meta={'Surface Brightess': 'table info'})
    

        
    def plot(self,color=None,ls='solid',smamax=np.inf):
        from photutils import EllipticalAperture
        # =============================================================================
        # Leemos las aperturas del fit
        # =============================================================================
        x0 = self.isolist.x0
        y0 = self.isolist.y0
        sma = self.isolist.sma
        smna = sma*(1 - self.isolist.eps)
        pa = self.isolist.pa
        # =============================================================================
        # Paleta de colores
        # =============================================================================
        if color==None:
            from matplotlib.pyplot import cm
            paleta=cm.rainbow(np.linspace(0,1,len(x0-1)))
        else:
            paleta=[color]*(len(x0-1))

        for i in range(1,len(x0)-1):
            if sma[i] <= smamax:
                ap = EllipticalAperture((x0[i], y0[i]), sma[i], smna[i], pa[i])
                ap.plot(color=paleta[i],ls=ls)





"""# My method for errors:

# error Poissoniano del bin - The error of the mean from ellipse.fit_image == (rms / sqrt(# data points))
error_poiss=isolist.int_err  

# error del cielo del bin
rms_sky=1. # Inventado, habría que calcularlo del cielo
error_sky=rms_sky/np.sqrt(isolist.ndata) 

# Total error
err_tot = np.sqrt(np.square(error_poiss)+np.square(error_sky))

# errorbars
up=isolist.intens+err_tot
down=isolist.intens-err_tot

# Pasamos los errores a magnitudes
error_up=[-2.5*(np.log10(up/isolist.intens))] 
error_down=[-2.5*(np.log10(isolist.intens/down))]
"""







# Clase vieja basada en Apertures
# class EllipticalProfile:
#     """ Descripcion
#     """
#     def __init__(self, image, Xc, Yc, nbins, a_min,npix, zp, rms_sky, ell, absorption, pixelscale, theta=0):
        
        
#         self.image=image
#         self.Xc = Xc
#         self.Yc = Yc
#         self.nbins = nbins
#         self.npix = npix
#         self.zp = zp
#         self.rms_sky = rms_sky
#         self.ell = ell
#         self.absorption = absorption
#         self.pixelscale = pixelscale
#         self.theta = theta
#         self.a_min = a_min 
        
#         #=======================================================================================================
#         #  Empezamos definiendo los parametros para las aperturas
#         #=======================================================================================================
#         import numpy as np
#         i=np.log10(a_min)
#         l=np.log10(self.npix)
        
#         step=np.logspace(i,l,self.nbins)
        
#         self.a=step
#         self.b=(1-self.ell)*self.a

#         #==============================================================================
#         # Aperturas
#         #==============================================================================
#         from photutils import EllipticalAperture
#         from photutils import EllipticalAnnulus
        
#         # La primera apertura es una ellipse y el resto son anillos concentricos
#         self.ap=[EllipticalAperture((self.Xc,self.Yc), a=self.a[0], b=self.b[0],theta=self.theta)]
        
#         for i in range(len(self.a)-1):
#             self.ap.append(EllipticalAnnulus((self.Xc,self.Yc), a_in=self.a[i],a_out=self.a[i+1], b_out=self.b[i+1],theta=self.theta))

#         #=========================================================================================================

#     def plot(self,image=None,color=None,ls='solid',cmap='viridis',vmin=None,vmax=None):
#         import matplotlib.pyplot as plt
#         # Si no le damos la imagen, pinta la imagen original del analisis.
#         if  np.shape(image)==() :
#             image=self.image

#         # Pintamos la imagen
#         plt.imshow(image,cmap=cmap,vmin=vmin,vmax=vmax,origin='lower')
        

#         # Paleta de colores
#         if color==None:
#             from matplotlib.pyplot import cm
#             paleta=cm.rainbow(np.linspace(0,1,len(self.a)))
#         else:
#             paleta=[color]*(len(self.x))
            
#         for i in range(len(self.ap)):
#             self.ap[i].plot(color=paleta[i],ls=ls)
        
#         #=========================================================================================================xw
            
#     def surfbrigth(self, method,subpixels=None,sigma=3):
#         from astropy.table import Table
#         from astropy.stats import sigma_clipped_stats
        
#         #==================================================================================
#         #  El primer elemento es la elipse central. El resto de las aperturas son anillos
#         #==================================================================================
        
#         sb = []
#         error_up = []
#         error_down = []

#         for i in range(len(self.ap)):
#             mask=self.ap[i].to_mask(method=method,subpixels=subpixels)
#             mask_data=mask.data
#             sec=mask.cutout(self.image)
#             sec_weight=sec*mask_data
#             sec_data=sec_weight[mask_data>0]            
            
#             """
#             plt.close('all')
#             plt.imshow(mask_data,vmin=0,vmax=10,origin='lower')
#             plt.show()
#             input()
#             """
#             mean= sigma_clipped_stats(sec_data, sigma=sigma)[0]
#             sb.append(-2.5*np.log10(mean) + self.zp + 5*np.log10(self.pixelscale) - self.absorption)
            
#             # error Poissoniano
#             error_poiss=np.sqrt(np.nanmean(np.square(mean - sec_data))/np.size(sec_data))      
#             # error bin
#             error_bin=self.rms_sky/np.sqrt(float(np.size(sec_data)))
        
#             # error total
#             up=mean+np.sqrt(np.square(error_poiss)+np.square(error_bin))
#             down=mean-np.sqrt(np.square(error_poiss)+np.square(error_bin))
#             # Pasamos los errores a magnitudes
#             error_up.append(-2.5*(np.log10(mean/up))) 
#             error_down.append(-2.5*(np.log10(down/mean)))
            
#             # Al finalizar definimos "r" como el punto medio de la anchura de la apertura
            
#             r=np.insert(self.a[1:]+np.diff(self.a)/2,0,self.a[0]/2)
#         #==================================================================================
#             vec = [r, sb, error_up, error_down]
#             names = ['r', 'sb','err up','err down']
#         return Table(vec, names=names, meta={'Surface Brighnes': 'table info'})
            
            
              
            
            
# =============================================================================
# Funcion para borrar secciones de las mascaras
# =============================================================================


def DeleteMaskRegion(mask, image, vmin= 0, vmax=0.3):
    """ Definicion
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from tkinter import messagebox #tkinter.TkVersion 8.6 
    from astropy.visualization import LogStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    
    plt.ion(),plt.show()

    # Establecemos los limites para asegura que siempre se pinten los mismos colores a las mismas regiones enmascaradas
    vmax_m=np.max(mask)

    remove=True
    plt.figure(1,figsize=(7,7))
    
    while remove==True:
        happy=False
        while happy == False:
            plt.cla()
            plt.imshow(mask,origin='lower',cmap='tab20c', aspect='equal',vmin=1,vmax=vmax_m)
            plt.tight_layout()
            plt.show()
        
            mask_new=np.copy(mask)
            image_new=np.copy(image)
            
            # Coger puntos de manera interactiva
            point=plt.ginput(n=1,timeout=0,show_clicks=False)
            coords = np.int_(point[0])
            value=mask[coords[1]][coords[0]]
            mask_new[mask_new==value]=0
            image_new[mask_new != 0]=np.nan

            
            plt.cla()
            #plt.imshow(mask_new,origin='lower',aspect='auto',vmin=vmin,vmax=vmax)
            
            # Create normalizer object
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

def HandMask(image, mask_value = True, vmin= 0, vmax=0.2):
    """ Definicion
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from tkinter import messagebox #tkinter.TkVersion 8.6 
    from astropy.visualization import LogStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    from photutils import CircularAperture

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
            
            # Coger puntos de manera interactiva preguntando cuantas mascaras se quieren (defecto 1)
            
            n=input('How many maks? (default 1): ')
            if n.isnumeric() == True:
                n=int(n)
            else:
                n=1
                
            radii=input('maks radii? (default 10): ')
            if radii.isnumeric() == True:
                radii=int(radii)
            else:
                radii=10
            
            # Buscamos las coordenadas e manera interactiva
            coords = plt.ginput(n=n,timeout=0,show_clicks=True)
            
            # Creamos una apertura circular en cada posicion
            ap = CircularAperture(coords, radii)
            
            shape = np.shape(image)
            
            # Hacemos un for sobre cada una de las posiciones
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
# Funcion para enmascarar una seccion eiptica en una posicion dada con el cursor
# =============================================================================

def EllipticalHandMask(image, mask_value = True, vmin= 0, vmax=0.2):
    """ Definicion
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
            
            # Coger puntos de manera interactiva preguntando cuantas mascaras se quieren (defecto 1)
            
            n=input('How many maks? (default 1): ')
            if n.isnumeric() == True:
                n=int(n)
            else:
                n=1
                
            a=input('semi-major axi ? (default 10): ')
            if a.isnumeric() == True:
                a=int(a)
            else:
                a=10

            b=input('semi-minor axi ? (default 5): ')
            if b.isnumeric() == True:
                b=int(b)
            else:
                b=5

            # Buscamos las coordenadas e manera interactiva
            coords = plt.ginput(n=n,timeout=0,show_clicks=True)
            
            # Creamos una apertura circular en cada posicion
            ap = EllipticalAperture(coords, a=a, b=b)
            
            shape = np.shape(image)
            
            # Hacemos un for sobre cada una de las posiciones
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

def Patchs(imagen, size):
    """ Definicion
    """
    y_l,x_l = np.shape(imagen)
    # Vamos a crear secciones de la imagen 
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
# LogNorm for imshow in og scale
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
# Look for GAMA groups within a RA, DEC region
# =============================================================================

def TableGAMAreg(ramin, ramax, decmin, decmax, tablepath):
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











#*****************************************************************************
#                             SB PROFILES FUNCTIONS 
#-----------------------------------------------------------------------------
#*****************************************************************************

# Package imports
import numpy as np
import matplotlib.pyplot as plt 

# Photutils modules
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import EllipticalAperture

# Astropy modules
from astropy.stats import sigma_clipped_stats        
from astropy.table import Table   
from astropy.io import fits
from astropy.wcs import WCS

# Scipy module
from scipy.optimize import minimize_scalar


# Versions used:
# numpy== 1.19.2
# matplotlib==3.3.2
# astropy==4.1
# photutils==1.0.1
# tkinter == 8.6

# To install specific versions, you can use:
# conda install astropy==4.1
# conda install photutils==1.0.1


#==============================================================================
# Ignore Warnings
#==============================================================================
import warnings
warnings.filterwarnings("ignore")

# Avoid warning if the data is masked, i.e., a single image with data + NaNs
warnings.warn("Input data contains invalid values (NaNs or infs), which were automatically clipped.") 


# =============================================================================
# Ancillary Functions
# =============================================================================

class Counts2Sb:
    def __init__(self, zp, pixscale, A=0, dimming=0, k_corr=0):
        '''
        Class to prepare the conversion from counts to surface brightness [mag/arcsec^2]

        Parameters
        ----------
        zp : float
            Zero point [mag]
        A : float, optional
            Galactic absorption or extinction for a given band. Default is 0.
        dimming : float, optional
            Dimming: 10*np.log10(1+z) [mag/arcsec^2]. 
            Because the dimming in linear scale is (1+z)**-4 [Calvi+2014, doi: 10.1088/0004-637X/796/2/102]. 
            Default is 0.
        k_corr : float, optional
            K correction [J. Loveday et al. (2012)]. For redshifted galaxies, the spectral shift changes the wavelength values as lambda*(1+z) along with the wavelength range.
            Default is 0.

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
        Function to convert counts to surface brightness [mag/arcsec^2] using the parameters from Counts2Sb.

        Parameters
        ----------
        counts : float or numpy array
            Values in counts.

        Returns
        -------
        Float or numpy array in [mag/arcsec^2].
        '''
        self.counts = counts
        # Uncomment the following line for debugging the conversion formula:
        # print(counts, -2.5*np.log10(counts) + self.zp + 5*np.log10(self.pixscale) - self.A - self.dimming - self.k_corr)
        return -2.5*np.log10(counts) + self.zp + 5*np.log10(self.pixscale) - self.A - self.dimming - self.k_corr 
    
    def do_errors(self, counts, counts_err):
        '''
        Function to convert count errors into surface brightness errors.

        Parameters
        ----------
        counts : float or numpy array
            Counts value.
        counts_err : float or numpy array
            Error in counts.

        Returns
        -------
        List containing the upper and lower errors in surface brightness.
        '''
        self.counts = counts 
        up = self.counts + counts_err
        down = self.counts - counts_err
        
        # Convert to magnitudes
        error_down = -2.5*np.log10(self.counts/up) 
        error_up = -2.5*np.log10(down/self.counts)
        return [error_up, error_down]

def sb2counts(sb, zp, pixscale, A=0, dimming=0, k_corr=0):
    '''
    Function to convert surface brightness [mag/arcsec^2] to counts.
    
    Parameters
    ----------
    sb : float or numpy array
        Surface brightness values [mag/arcsec^2].
    
    Returns
    -------
    Float or numpy array representing counts.
    '''
    return 10**(-0.4*(sb - zp - 5.*np.log10(pixscale) + A + dimming + k_corr) )
    
    
def SecImData(ap, image, method, subpixels=None):
    """
    Extracts a section of a 2D array using any astropy.aperture shape.

    Parameters
    ----------
    ap : astropy.aperture class
        Aperture object defining the geometry.
    image : numpy array
        Data from which to extract the section.
    method : str
         Interpolation method to use:
         'exact': calculates the exact intersection of the aperture with each pixel.
         'center': a pixel is considered entirely inside or outside the aperture based on its center.
         'subpixel': divides pixels into subpixels to decide based on their centers (requires specifying the number of subpixels).
         (Note: 'center' and 'subpixel' methods are faster but less precise.)
    subpixels : int, optional
        Number of subpixels to use with the 'subpixel' method. Default is None.

    Returns
    -------
    sec_data : 1D numpy array
        1D array of data within the aperture, weighted as specified by the method,
        with NaN values removed.
    """
    mask = ap.to_mask(method=method, subpixels=subpixels)
    mask_data = mask.data
    sec = mask.cutout(image)
    sec_weight = sec * mask_data
    sec_data = sec_weight[mask_data > 0]
    # Remove NaN values
    sec_data = sec_data[np.isfinite(sec_data)]
    
    return sec_data


def com(data, Xcoords, Ycoords):
    """
    Function to compute the center of light (or mass) for a set of positions in 
    an astronomical image (or any array). For example, it can determine the center
    of light of a system of three galaxies given their central coordinates.

    Parameters
    ----------
    data : np.array
        Image containing the sources.
    Xcoords : list
        X coordinates of the sources (as used in ds9).
    Ycoords : list 
        Y coordinates of the sources (as used in ds9).

    Returns
    -------
    list
        [x, y] coordinates of the center of light (or mass) in Python convention.
    """
    weights_list = []

    for ii in np.arange(len(Xcoords)):
        weights_list.append(data[Ycoords[ii], Xcoords[ii]])

    weights = np.array(weights_list)
    
    # Calculate the center of light using a weighted average
    com_x = sum(Xcoords * weights) / sum(weights)
    com_y = sum(Ycoords * weights) / sum(weights)

    return [com_x, com_y]



class Center:
    # Adjustment to find the center
    def __init__(self, image, Xc, Yc, d):
        self.image = image
        self.Xc = Xc
        self.Yc = Yc
        self.d = d

        xmas = np.int(Xc + d)
        xmenos = np.int(Xc - d)  
        ymas = np.int(Yc + d)
        ymenos = np.int(Yc - d)
        
        # Build a section to calculate the centroid
        self.sec = image[ymenos:ymas, xmenos:xmas]
        
        from photutils import centroid_2dg  # Also available: centroid_1dg, centroid_com, etc.
        
        # The positions of the fitted centroids
        x_s, y_s = centroid_2dg(self.sec)
        
        self.x = x_s + self.Xc - d
        self.y = y_s + self.Yc - d

    def find(self):              
        return (self.x, self.y)
    
    def plot(self):
        p_x = np.sum(self.sec, 0)
        p_y = np.sum(self.sec, 1)
        ly = np.arange(len(p_y)) + self.Yc - self.d
        lx = np.arange(len(p_x)) + self.Xc - self.d
        
        plt.close('all')
        plt.figure(1, figsize=(8,4))

        # Profile along the X axis
        plt.subplot(121)
        plt.title('X direction')
        plt.step(lx, p_x)
        plt.axvline(self.x, color='r', ls='dashed')
        plt.xlabel('pix')
        plt.ylabel('Counts')
        
        # Profile along the Y axis
        plt.subplot(122)
        plt.title('Y direction')
        plt.step(ly, p_y)
        plt.axvline(self.y, color='r', ls='dashed')
        plt.xlabel('pix')
        plt.yticks([], [])
        
        plt.show()


#********************************************************************************
# RADIAL CENTRAL PROFILE 
#********************************************************************************

#==============================================================================
# Radial profiles. Root function
#==============================================================================

class CentralProfile:
    def __init__(self, image, Xc, Yc, nbins, npix, height, zp, pixscale, A=0, rms_sky=0, orientation='horizontal', delta=0):
        """
        Main class definition to extract a surface brightness profile using rectangular apertures over an image.
        The apertures are built here.
        
        Parameters
        ----------
        image : 2D numpy array
            Data from which to extract the profile.
        Xc : float
            X center coordinate in pixels.
        Yc : float
            Y center coordinate in pixels.
        nbins : int
            Number of bins in the profile.
        npix : int
            The total width of the profile in pixels; the number of pixels up to which the profile is extracted.
        height : float
            The full height of the rectangular apertures in pixels.
        zp : float
            Zero point [mag/arcsec^2]. Default is 0.
        A : float
            Galactic absorption or extinction for a given band. Default is 0.
        rms_sky : float
            Sky RMS to include as a source of error. Default is 0.
        orientation : str, optional
            Orientation of the profile. Default is 'horizontal'. 'horizontal' is along the galaxy's major axis; 'vertical' along the minor axis.
        delta : float, optional
            Pixel offset from the galaxy center at which to extract the profile. Default is 0, i.e., along the main axis.

        Returns
        -------
        Instance of the class with the apertures already defined.
        """
        # Geometric parameters
        self.image = image
        self.Xc = Xc
        self.Yc = Yc
        self.nbins = nbins
        self.npix = npix
        self.height = height
        self.delta = delta
        self.orientation = orientation
        
        # Data parameters
        self.zp = zp
        self.A = A
        self.rms_sky = rms_sky
        self.pixscale = pixscale
        
        #=======================================================================================
        # Apertures definition
        #=======================================================================================
        import numpy as np
        l = np.log10(self.npix)
        step = np.logspace(0, l, self.nbins)

        # Dictionary for all parameters depending on whether the aperture is vertical or horizontal.
        self.param = {
            'horizontal': {
                'w': np.diff(step),
                'h': np.ones(len(np.diff(step)))*self.height,
                'x': step[1:] - (np.diff(step))/2 - step[0] - np.diff(step)[0]/2,
                'y': np.zeros(len(np.diff(step))) + self.delta,
                'signo': (-1, 1)
            },
            'vertical': {
                'w': np.ones(len(np.diff(step)))*self.height,
                'h': np.diff(step),
                'x': np.zeros(len(np.diff(step))) + self.delta,
                'y': step[1:] - (np.diff(step))/2 - step[0] - np.diff(step)[0]/2,
                'signo': (1, -1)
            }
        }

        w = self.param[orientation]['w']
        h = self.param[orientation]['h']
        x = self.param[orientation]['x']
        y = self.param[orientation]['y']
        self.signo = self.param[orientation]['signo']

        #==============================================================================
        # Build apertures
        #==============================================================================
        from photutils import RectangularAperture

        self.ap_c = RectangularAperture([x[0] + self.Xc, y[0] + self.Yc], w=w[0], h=h[0], theta=0)  # Central aperture
        self.ap_r = []  # Apertures to the right
        self.ap_l = []  # Apertures to the left
        
        for x0, y0, w0, h0 in zip(x[1:], y[1:], w[1:], h[1:]):
            self.ap_r.append(RectangularAperture((x0 + self.Xc, y0 + self.Yc), w=w0, h=h0, theta=0))
            self.ap_l.append(RectangularAperture((self.signo[0]*x0 + self.Xc, self.signo[1]*y0 + self.Yc), w=w0, h=h0, theta=0))

        self.w = w
        
        #==========================================================================================
        
    def plot(self, color=None, ls='solid', alpha=1, lw=1, n_max=np.inf):
        """
        Plot the profile apertures.

        Parameters
        ----------
        color : str, optional
            Colour of the apertures in the plot. Default is None.
        ls : str, optional
            Linestyle (e.g., "solid", "dotted", "dashed", "dashdot", etc.). Default is 'solid'.
        alpha : float, optional
            Transparency of the plotted lines. Default is 1.
        lw : float, optional
            Line width. Default is 1.
        n_max : int, optional
            Maximum number of apertures to plot. Default is infinity.

        Returns
        -------
        Displays a plot of the apertures.
        """
                
        # Colour palette
        if color is None:
            from matplotlib.pyplot import cm
            paleta = cm.rainbow(np.linspace(0, 1, len(self.ap_r) + 1))
        else:
            paleta = [color] * (len(self.ap_r) + 1)
        
        # Plot the central aperture
        self.ap_c.plot(color=paleta[0], ls=ls, alpha=alpha, lw=lw)
        x, y = self.ap_c.positions
        # plt.text(x, y, 0, ha='center', va='bottom')
           
        for i, (ap_r0, ap_l0) in enumerate(zip(self.ap_r, self.ap_l)):
            if i < n_max:
                # Right aperture
                ap_r0.plot(color=paleta[i+1], ls=ls, alpha=alpha, lw=lw)
                x, y = ap_r0.positions
                # plt.text(x, y, i+1, ha='center', va='bottom', color=paleta[i+1])
                # Left aperture
                ap_l0.plot(color=paleta[i+1], ls=ls, alpha=alpha, lw=lw)
                x, y = ap_l0.positions
                # plt.text(x, y, i+1, ha='center', va='bottom', color=paleta[i+1])
                
    def saveDS9(self, fname):
        """
        Write the apertures to a file as DS9 region definitions.

        Parameters
        ----------
        fname : str
            Filename (without extension).

        Returns
        -------
        A .reg file saved in the current folder.
        """
        # Header for all regions
        l1 = '# Region file format: DS9 version 4.1'
        l2 = 'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
        l3 = 'image'
        
        file = open(fname + '.reg', 'w')
        file.write(l1 + "\n")
        file.write(l2 + "\n")
        file.write(l3 + "\n")
        
        # Write the line for the central aperture
        x, y = self.ap_c.positions
        w = self.ap_c.w
        h = self.ap_c.h
        
        l = 'box(' + str(x) + ',' + str(y) + ',' + str(w) + ',' + str(h) + ') # text={0} dash=1'
        file.write(l + "\n")

        # Now write the lines for the remaining apertures
        for i, (ap_r0, ap_l0) in enumerate(zip(self.ap_r, self.ap_l)):          
            # Right apertures
            x, y = ap_r0.positions
            w = ap_r0.w
            h = ap_r0.h
        
            l_r = 'box(' + str(x) + ',' + str(y) + ',' + str(w) + ',' + str(h) + ') # text={' + str(i+1) + '} dash=1'
            file.write(l_r + "\n")
            
            # Left apertures
            x, y = ap_l0.positions
            w = ap_l0.w
            h = ap_l0.h
        
            l_l = 'box(' + str(x) + ',' + str(y) + ',' + str(w) + ',' + str(h) + ') # text={' + str(i+1) + '} dash=1'
            file.write(l_l + "\n")
    
        file.close()
    
    
    def surfbrigth(self, method, subpixels=None, sigma=3):
        """
        Extracts the surface brightness profile from the apertures defined in the CentralProfile class.

        Parameters
        ----------
        method : str
            Interpolation method:
            'exact': calculates the exact intersection of the aperture with each pixel.
            'center': a pixel is considered entirely inside or outside the aperture depending on whether its center is inside.
            'subpixel': pixels are divided into subpixels and classified based on their centers. (Note: 'center' and 'subpixel' are faster, but less precise.)
        subpixels : int, optional
            Number of subpixels for the 'subpixel' method. Default is None.
        sigma : float, optional
            Sigma value for the sigma-clipping rejection algorithm in the statistics. Default is 3.

        Returns
        -------
        astropy.table.Table
            Table with profile values in counts and surface brightness, and their errors.
        """
        
        from astropy.table import Table
        from astropy.stats import sigma_clipped_stats  
        from astropy.stats import sigma_clip
        # Auxiliary class that defines the conversion from counts to surface brightness
        c2sb = Counts2Sb(self.zp, self.pixscale, self.A)

        # =============================================================================
        # The first element for both the right and left sides is the central aperture.
        # =============================================================================
        # ---------------------------------------------------------------------------
        # Spatial variables initialization
        # ---------------------------------------------------------------------------
        x = [0]
        y = [0]
        
        # Extract data from the central aperture section
        sec_data = SecImData(self.ap_c, self.image, method, subpixels=subpixels)
        # Compute the sigma-clipped median of the central aperture
        median = sigma_clipped_stats(sec_data, sigma=sigma)[1]  
        
        # =============================================================================
        # Surface brightness calculation
        # =============================================================================
        # Counts for central aperture
        c_r = [median]
        c_l = [median]
        c = [median]
        
        # Sum of flux in the bin for the central aperture
        sum_ct = np.nansum(sec_data)
        sum_sb = [c2sb.do(sum_ct)]
        
        # Convert counts to magnitudes (surface brightness)
        sb_r = [c2sb.do(median)]
        sb_l = [c2sb.do(median)]
        sb = [c2sb.do(median)]
        
        # =============================================================================
        # Errors Calculation
        # =============================================================================        
        # Poisson error
        error_poiss = np.sqrt(np.nanmean(np.square(median - sec_data)) / np.size(sec_data))      
        # Error from sky RMS (bin error)
        error_bin = self.rms_sky / np.sqrt(float(np.size(sec_data)))
    
        # Total error in counts
        err_t = np.sqrt(np.square(error_poiss) + np.square(error_bin))
        c_err = [err_t]
        c_err_r = [err_t]
        c_err_l = [err_t]
 
        # For magnitude conversion, compute upper and lower counts
        up = median + err_t
        down = median - err_t
        
        # Convert counts errors to surface brightness errors
        error_down_r = [-2.5 * np.log10(median / up)]
        error_up_r = [-2.5 * np.log10(down / median)]

        error_down_l = [-2.5 * np.log10(median / up)]
        error_up_l = [-2.5 * np.log10(down / median)]
        
        error_down = [-2.5 * np.log10(median / up)]
        error_up = [-2.5 * np.log10(down / median)]
        
        # =============================================================================
        # Bin size in pixels for the central aperture       
        # =============================================================================
        bin_size = [float(np.size(sec_data))]
        
        # =============================================================================
        # Process the remaining apertures (non-central)
        # =============================================================================
        for ap_r0, ap_l0 in zip(self.ap_r, self.ap_l):
            # Compute spatial offset relative to the central aperture
            x.append(ap_r0.positions[0] - self.ap_c.positions[0])
            y.append(ap_r0.positions[1] - self.ap_c.positions[1])

            # ---------------------------------------------------------------------------
            # Right aperture processing
            # ---------------------------------------------------------------------------
            sec_data_r = SecImData(ap_r0, self.image, method, subpixels=subpixels) 
            # Sigma-clipped median for the right aperture
            median_r = sigma_clipped_stats(sec_data_r, sigma=sigma)[1]

            # ---------------------------------------------------------------------------
            # Left aperture processing
            # ---------------------------------------------------------------------------
            sec_data_l = SecImData(ap_l0, self.image, method, subpixels=subpixels)
            # Sigma-clipped median for the left aperture
            median_l = sigma_clipped_stats(sec_data_l, sigma=sigma)[1]
            
            # Average counts from both sides
            mean_rl = np.nanmean([median_r, median_l])

            # Update counts for right, left and average
            c_r.append(median_r)
            c_l.append(median_l)
            c.append(mean_rl)
            
            # Convert counts to surface brightness (magnitudes) for each side and average
            sb_r.append(c2sb.do(median_r))
            sb_l.append(c2sb.do(median_l))
            sb.append(c2sb.do(mean_rl))
        
            # Sum flux for the bin in each side
            sum_ct_r = np.nansum(sec_data_r)
            sum_sb_r = c2sb.do(sum_ct_r)
            sum_ct_l = np.nansum(sec_data_l)
            sum_sb_l = c2sb.do(sum_ct_l)
            # Average the summed flux from both sides
            sum_ct = np.nanmean([sum_ct_r, sum_ct_l])
            sum_sb.append(c2sb.do(sum_ct)) 
            
            # =============================================================================
            # Error estimation for each aperture
            # =============================================================================
            # Poisson error for right and left apertures
            error_poiss_r = np.sqrt(np.nanmean(np.square(median_r - sec_data_r)) / np.size(sec_data_r)) 
            error_poiss_l = np.sqrt(np.nanmean(np.square(median_l - sec_data_l)) / np.size(sec_data_l))
            error_poiss = np.sqrt(np.square(error_poiss_r) + np.square(error_poiss_l))
            
            # Bin errors for right and left apertures
            error_bin_r = self.rms_sky / np.sqrt(float(np.size(sec_data_r)))
            error_bin_l = self.rms_sky / np.sqrt(float(np.size(sec_data_l)))
            error_bin = np.sqrt(np.square(error_bin_r) + np.square(error_bin_l))

            # Total error in counts for each side
            err_r = np.sqrt(np.square(error_poiss_r) + np.square(error_bin_r))
            c_err_r.append(err_r)
            
            err_l = np.sqrt(np.square(error_poiss_l) + np.square(error_bin_l))
            c_err_l.append(err_l)
            
            # Total error for the combined aperture
            err_t = np.sqrt(np.square(error_poiss) + np.square(error_bin))
            c_err.append(err_t)

            # Convert counts error to surface brightness error for each side and overall
            up_r = median_r + err_r
            down_r = median_r - err_r
            
            up_l = median_l + err_l
            down_l = median_l - err_l
            
            up = mean_rl + err_t
            down = mean_rl - err_t
            
            error_down_r.append(-2.5 * np.log10(median_r / up_r)) 
            error_up_r.append(-2.5 * np.log10(down_r / median_r))

            error_down_l.append(-2.5 * np.log10(median_l / up_l)) 
            error_up_l.append(-2.5 * np.log10(down_l / median_l))

            error_down.append(-2.5 * np.log10(mean_rl / up)) 
            error_up.append(-2.5 * np.log10(down / mean_rl))
                        
            # =============================================================================
            # Bin size in pixels for each non-central aperture         
            # =============================================================================
            bin_size.append(float(np.size(sec_data_r)))     
            
        # =============================================================================
        # Replace any NaN values with 0
        # =============================================================================
        error_up = np.nan_to_num(error_up)
        error_down = np.nan_to_num(error_down)

        error_down_r = np.nan_to_num(error_down_r)
        error_up_r = np.nan_to_num(error_up_r)

        error_down_l = np.nan_to_num(error_down_l)
        error_up_l = np.nan_to_num(error_up_l)

        # =============================================================================
        # Define the spatial variable dictionary and convert to arcseconds
        # =============================================================================
        self.param['horizontal'] = {'spatial': ('r', x)}
        self.param['vertical'] = {'spatial': ('z', y)}

        spatial_name = self.param[self.orientation]['spatial'][0]
        spatial = np.array(self.param[self.orientation]['spatial'][1]) * self.pixscale
        
        bins_width = self.w * self.pixscale
        
        # Create the table data vector and corresponding column names
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

        return Table(vec, names=names, meta={'Surface Brighnes': 'table info, u.arcsec'})
        
        
#==============================================================================
# Shifted radial profiles
#==============================================================================

class ShiftedProfile:
    
    def __init__(self, image, Xc, Yc, nbins, npix, height, zp, pixscale, delta, rms_sky=0, A=0, orientation='horizontal'):
        '''
        Class definition to extract a surface brightness profile with an offset from the main axis of the galaxy
        using rectangular apertures over an image. The apertures are built here.
        
        Parameters
        ----------
        image : 2D numpy array
            Data from which to extract the profile.
        Xc : float
            X center coordinate in pixels.
        Yc : float
            Y center coordinate in pixels.
        nbins : int
            Number of bins in the profile.
        npix : int
            Total width of the profile in pixels (number of pixels up to which the profile is extracted).
        height : float
            Height of the rectangular apertures in pixels.
        zp : float
            Zero point [mag/arcsec^2].
        delta : float
            Offset in pixels from the galaxy center at which to extract the profile.
        rms_sky : float
            Sky RMS to include as a source of error.
        A : float
            Galactic absorption or extinction for a given band.
        orientation : str, optional
            Orientation of the profile. 'horizontal' for the galaxy's major axis; 'vertical' for the minor axis. Default is 'horizontal'.

        Returns
        -------
        None.
        '''
        self.image = image
        self.Xc = Xc
        self.Yc = Yc
        self.nbins = nbins
        self.npix = npix
        self.height = height
        self.zp = zp
        self.A = A
        self.rms_sky = rms_sky
        self.orientation = orientation
        self.delta = delta
        self.pixscale = pixscale
        # Instantiate the main CentralProfile class for the upward offset
        self.up = CentralProfile(image, Xc, Yc, nbins, npix, height, zp, pixscale, A, rms_sky, orientation, delta)
        # Instantiate the main CentralProfile class for the downward offset (negative delta)
        self.down = CentralProfile(image, Xc, Yc, nbins, npix, height, zp, pixscale, A, rms_sky, orientation, -delta)

    def plot(self, color=None, ls='solid', alpha=1, lw=1, n_max=np.inf):
        '''
        Plot the profile apertures for both the upward and downward profiles.

        Parameters
        ----------
        color : str, optional
            Colour for the apertures in the plot. Default is None.
        ls : str, optional
            Linestyle (e.g., "solid", "dotted", "dashed", "dashdot", etc.). Default is 'solid'.

        Returns
        -------
        Displays a plot of the apertures.
        '''
        self.up.plot(color=color, ls=ls, alpha=alpha, lw=lw, n_max=n_max)
        self.down.plot(color=color, ls=ls, alpha=alpha, lw=lw, n_max=n_max)
    
    def saveDS9(self, fname):
        '''
        Write the apertures to a file as DS9 region definitions.

        Parameters
        ----------
        fname : str
            Filename (without extension).

        Returns
        -------
        A .reg file saved in the current folder.
        '''
        self.up.saveDS9(fname + '_up')
        self.down.saveDS9(fname + '_down')        
        

    def surfbrigth(self, method, subpixels=None, sigma=3):
        '''
        Extracts the surface brightness profile from the apertures defined in the ShiftedProfile class.

        Parameters
        ----------
        method : str
            Interpolation method as described previously.
        subpixels : int, optional
            Number of subpixels for the 'subpixel' method. Default is None.
        sigma : float, optional
            Sigma for the sigma-clipping algorithm. Default is 3.

        Returns
        -------
        astropy.table.Table
            Table with profile values in counts and surface brightness, and their errors.
        '''
        
        # Auxiliary class for converting counts to surface brightness
        c2sb = Counts2Sb(self.zp, self.pixscale, self.A)

        # Obtain the profile table from the upward and downward profiles
        t_1 = self.up.surfbrigth(method, subpixels, sigma)
        t_2 = self.down.surfbrigth(method, subpixels, sigma)
        
        # Retrieve the spatial coordinate name and values from one of the tables (they are the same)
        spatial = t_1.columns[0]
        spatial_name = t_1.colnames[0]
       
        # Upward/downward apertures:
        # Counts
        ct_up = t_1['ct']
        ct_down = t_2['ct']
        # Surface brightness for up/down apertures
        sb_up = t_1['sb']
        sb_down = t_2['sb']
        
        # Right apertures:
        # Counts
        ct_up_r = t_1['ct right']
        ct_down_r = t_2['ct right']
        # Average counts for right apertures
        ct_r = np.nanmean([ct_up_r, ct_down_r], 0)
        # Surface brightness for right apertures
        sb_r = c2sb.do(ct_r)
        sb_up_r = t_1['sb right']
        sb_down_r = t_2['sb right']
        
        # Left apertures:
        # Counts
        ct_up_l = t_1['ct left']
        ct_down_l = t_2['ct left']
        # Average counts for left apertures
        ct_l = np.nanmean([ct_up_l, ct_down_l], 0)
        # Surface brightness for left apertures
        sb_l = c2sb.do(ct_l)
        sb_up_l = t_1['sb left']
        sb_down_l = t_2['sb left']
        
        # Average over all apertures:
        # Counts
        ct_mean = np.nanmean([ct_up_r, ct_down_r, ct_up_l, ct_down_l], 0)
        # Surface brightness for the average
        sb_mean = c2sb.do(ct_mean)

        # =============================================================================
        # Errors Calculation
        # =============================================================================
        # Errors in counts for up/down apertures
        ct_err_up = t_1['ct err']
        ct_err_down = t_2['ct err']        

        # Right apertures errors
        ct_err_up_r = t_1['ct err right']
        ct_err_down_r = t_2['ct err right']  
        ct_err_r = np.sqrt(np.square(ct_err_up_r) + np.square(ct_err_down_r))

        # Left apertures errors
        ct_err_up_l = t_1['ct err left']
        ct_err_down_l = t_2['ct err left']   
        ct_err_l = np.sqrt(np.square(ct_err_up_l) + np.square(ct_err_down_l))

        # Total counts error (up and down)
        ct_err = np.sqrt(np.square(ct_err_up) + np.square(ct_err_down))        
        
        #######################################################################
        # Convert counts errors to surface brightness errors:
        # For up/down apertures
        err_up_barDown = -2.5 * np.log10(ct_up / (ct_up + ct_err_up))
        err_up_barUp = -2.5 * np.log10((ct_up - ct_err_up) / ct_up)
        
        err_down_barDown = -2.5 * np.log10(ct_down / (ct_down + ct_err_down))
        err_down_barUp = -2.5 * np.log10((ct_down - ct_err_down) / ct_down)

        # For right apertures
        err_up_r_barDown = -2.5 * np.log10(ct_up_r / (ct_up_r + ct_err_up_r))
        err_up_r_barUp = -2.5 * np.log10((ct_up_r - ct_err_up_r) / ct_up_r)

        err_down_r_barDown = -2.5 * np.log10(ct_down_r / (ct_down_r + ct_err_down_r))
        err_down_r_barUp = -2.5 * np.log10((ct_down_r - ct_err_down_r) / ct_down_r)

        err_r_barDown = -2.5 * np.log10(ct_r / (ct_r + ct_err_r))
        err_r_barUp = -2.5 * np.log10((ct_r - ct_err_r) / ct_r)        

        # For left apertures
        err_up_l_barDown = -2.5 * np.log10(ct_up_l / (ct_up_l + ct_err_up_l))
        err_up_l_barUp = -2.5 * np.log10((ct_up_l - ct_err_up_l) / ct_up_l)

        err_down_l_barDown = -2.5 * np.log10(ct_down_l / (ct_down_l + ct_err_down_l))
        err_down_l_barUp = -2.5 * np.log10((ct_down_l - ct_err_down_l) / ct_down_l)

        err_l_barDown = -2.5 * np.log10(ct_l / (ct_l + ct_err_l))
        err_l_barUp = -2.5 * np.log10((ct_l - ct_err_l) / ct_l)          
        
        # Total surface brightness error for the combined aperture
        err_barDown = -2.5 * np.log10(ct_mean / (ct_mean + ct_err))
        err_barUp = -2.5 * np.log10((ct_mean - ct_err) / ct_mean)
        
        # =============================================================================
        # Create a table with all the data
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
# Projection onto a new image
# =============================================================================

class ProjectCentralProfile:
    def __init__(self, image, profile, wcs_ori, wcs, zp, pixscale, A=0, rms_sky=0, dx=0, dy=0):
        # Dictionary for all parameters depending on whether the aperture is vertical or horizontal.
        self.param = {'horizontal': {'signo': (-1, 1)},
                      'vertical': {'signo': (1, -1)}}

        # Geometric parameters
        self.orientation = profile.orientation
        
        # World Coordinate System (WCS)
        self.wcs_ori = wcs_ori 
        self.wcs = wcs
        
        # Data parameters
        self.image = image
        self.zp = zp
        self.A = A
        self.rms_sky = rms_sky
        self.pixscale = pixscale
        self.profile = profile

        # =============================================================================
        # Apertures definition
        # =============================================================================
        # Read the apertures from the original class and convert them to pixel coordinates
        # of the other image.
        self.ap_c = profile.ap_c.to_sky(wcs_ori).to_pixel(wcs)
        self.ap_c.positions[0] = self.ap_c.positions[0] - dx
        self.ap_c.positions[1] = self.ap_c.positions[1] - dy
        
        self.ap_r = []  # Apertures to the right
        self.ap_l = []  # Apertures to the left

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

    def plot(self, color=None, ls='solid'):
        """
        Plot the profile apertures.

        Parameters
        ----------
        color : str, optional
            Colour of the apertures in the plot. Default is None.
        ls : str, optional
            Linestyle ("solid", "dotted", "dashed", "dashdot", etc.). Default is 'solid'.
        cent_ap : bool, optional
            Plot a mark at the central position of the aperture. (Not used here.)
        
        Returns
        -------
        Displays a plot of the apertures.
        """
                
        # Colour palette
        if color is None:
            from matplotlib.pyplot import cm
            paleta = cm.rainbow(np.linspace(0, 1, len(self.ap_r) + 1))
        else:
            paleta = [color] * (len(self.ap_r) + 1)
        
        # Plot the central aperture
        self.ap_c.plot(color=paleta[0], ls=ls)
        x, y = self.ap_c.positions
        plt.text(x, y, 0, ha='center', va='bottom')
           
        for i, (ap_r0, ap_l0) in enumerate(zip(self.ap_r, self.ap_l)):
            # Right aperture
            ap_r0.plot(color=paleta[i + 1], ls=ls)
            x, y = ap_r0.positions
            plt.text(x, y, i + 1, ha='center', va='bottom', color=paleta[i + 1])
            # Left aperture
            ap_l0.plot(color=paleta[i + 1], ls=ls)
            x, y = ap_l0.positions
            plt.text(x, y, i + 1, ha='center', va='bottom', color=paleta[i + 1])
                
    def saveDS9(self, fname):
        """
        Writes the apertures to a file as DS9 region definitions.

        Parameters
        ----------
        fname : str
            File name.

        Returns
        -------
        A .reg file saved in the current folder.
        """
        # Header common to all regions
        l1 = '# Region file format: DS9 version 4.1'
        l2 = 'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
        l3 = 'image'
        
        file = open(fname + '.reg', 'w')
        file.write(l1 + "\n")
        file.write(l2 + "\n")
        file.write(l3 + "\n")
        
        # Write the line for the central aperture
        x, y = self.ap_c.positions
        w = self.ap_c.w
        h = self.ap_c.h
        
        l = 'box(' + str(x) + ',' + str(y) + ',' + str(w) + ',' + str(h) + ') # text={0} dash=1'
        file.write(l + "\n")

        # Now write for the rest of the apertures
        for i, (ap_r0, ap_l0) in enumerate(zip(self.ap_r, self.ap_l)):
            # Right apertures
            x, y = ap_r0.positions
            w = ap_r0.w
            h = ap_r0.h
        
            l_r = 'box(' + str(x) + ',' + str(y) + ',' + str(w) + ',' + str(h) + ') # text={' + str(i + 1) + '} dash=1'
            file.write(l_r + "\n")
            
            # Left apertures
            x, y = ap_l0.positions
            w = ap_l0.w
            h = ap_l0.h
        
            l_l = 'box(' + str(x) + ',' + str(y) + ',' + str(w) + ',' + str(h) + ') # text={' + str(i + 1) + '} dash=1'
            file.write(l_l + "\n")
    
        file.close()
    
    
    
    def surfbrigth(self, method, subpixels=None, sigma=3):
        """
        Extracts the surface brightness profile from the apertures defined in the CentralProfile class.

        Parameters
        ----------
        method : str
            Interpolation method:
            'exact': calculates the exact intersection of the aperture with each pixel.
            'center': a pixel is considered entirely inside or outside the aperture based on its centre.
            'subpixel': pixels are divided into subpixels and classified based on their centres. ('center' and 'subpixel' are faster, but less precise.)
        subpixels : int, optional
            Number of subpixels if 'subpixel' method is used. Default is None.
        sigma : float, optional
            Sigma value for the sigma-clipping rejection algorithm in the statistics. Default is 3.

        Returns
        -------
        astropy.table.Table
            Table with profile values in counts and surface brightness, and their errors.
        """
        
        from astropy.table import Table
        from astropy.stats import sigma_clipped_stats        
        # Auxiliary class that defines the conversion from counts to surface brightness
        c2sb = Counts2Sb(self.zp, self.pixscale, self.A)

        # =============================================================================
        # The first element for both the right and left sides is the central aperture.
        # =============================================================================
        # ---------------------------------------------------------------------------
        # Spatial variables initialization
        # ---------------------------------------------------------------------------
        x = [0]
        y = [0]
        
        # Extract section within the central aperture
        sec_data = SecImData(self.ap_c, self.image, method, subpixels=subpixels)
        # Compute sigma-clipped median for the central aperture
        median = sigma_clipped_stats(sec_data, sigma=sigma)[1]  
        
        
        # =============================================================================
        # Surface brightness calculation
        # =============================================================================
        # Counts
        c_r = [median]
        c_l = [median]
        c = [median]
        
        # Convert total flux in counts to surface brightness
        sb_r = [c2sb.do(median)]
        sb_l = [c2sb.do(median)]
        sb = [c2sb.do(median)]
        
        # =============================================================================
        # Errors
        # =============================================================================        
        # Poisson error
        error_poiss = np.sqrt(np.nanmean(np.square(median - sec_data)) / np.size(sec_data))      
        # Bin error based on sky RMS
        error_bin = self.rms_sky / np.sqrt(float(np.size(sec_data)))
    
        # Total error in counts
        err_t = np.sqrt(np.square(error_poiss) + np.square(error_bin))
        c_err = [err_t]
 
        # Compute upper and lower counts for magnitude conversion
        up = median + err_t
        down = median - err_t
        
        # Convert count errors to surface brightness errors
        error_down = [-2.5 * np.log10(median / up)]
        error_up = [-2.5 * np.log10(down / median)]
        
        # =============================================================================
        # Process the remaining apertures
        # =============================================================================
        for ap_r0, ap_l0 in zip(self.ap_r, self.ap_l):
            # Spatial variables: project the spatial coordinates to the original coordinate system.
            x0 = ap_r0.to_sky(self.wcs).to_pixel(self.wcs_ori).positions[0] - self.ap_c.to_sky(self.wcs).to_pixel(self.wcs_ori).positions[0]
            y0 = ap_r0.to_sky(self.wcs).to_pixel(self.wcs_ori).positions[1] - self.ap_c.to_sky(self.wcs).to_pixel(self.wcs_ori).positions[1]
            
            x.append(x0)
            y.append(y0)

            # Right apertures
            sec_data_r = SecImData(ap_r0, self.image, method, subpixels=subpixels) 
            # Compute sigma-clipped median for right aperture
            median_r = sigma_clipped_stats(sec_data_r, sigma=sigma)[1]

            # Left apertures
            sec_data_l = SecImData(ap_l0, self.image, method, subpixels=subpixels)
            # Compute sigma-clipped median for left aperture
            median_l = sigma_clipped_stats(sec_data_l, sigma=sigma)[1]
            
            # Average the medians from right and left apertures
            mean_rl = np.nanmean([median_r, median_l])

            # =============================================================================
            # Surface brightness for each aperture
            # =============================================================================
            # Counts
            c_r.append(median_r)
            c_l.append(median_l)
            c.append(mean_rl)
            
            # Magnitudes
            sb_r.append(c2sb.do(median_r))
            sb_l.append(c2sb.do(median_l))
            sb.append(c2sb.do(mean_rl))
        
            # =============================================================================
            # Errors
            # =============================================================================
            # Poisson error for right and left apertures
            error_poiss_r = np.sqrt(np.nanmean(np.square(median_r - sec_data_r)) / np.size(sec_data_r)) 
            error_poiss_l = np.sqrt(np.nanmean(np.square(median_l - sec_data_l)) / np.size(sec_data_l))
            error_poiss = np.sqrt(np.square(error_poiss_r) + np.square(error_poiss_l))
            
            # Bin errors for right and left apertures
            error_bin_r = self.rms_sky / np.sqrt(float(np.size(sec_data_r)))
            error_bin_l = self.rms_sky / np.sqrt(float(np.size(sec_data_l)))
            error_bin = np.sqrt(np.square(error_bin_r) + np.square(error_bin_l))

            # Total error in counts for this bin
            err_t = np.sqrt(np.square(error_poiss) + np.square(error_bin))
            c_err.append(err_t)

            # Magnitude errors: compute upper and lower values
            up = mean_rl + err_t
            down = mean_rl - err_t
            
            error_down.append(-2.5 * np.log10(mean_rl / up))
            error_up.append(-2.5 * np.log10(down / mean_rl))
            
        # =============================================================================
        # Replace any NaN values with 0
        # =============================================================================
        error_up = np.nan_to_num(error_up)
        error_down = np.nan_to_num(error_down)

        # =============================================================================
        # Define the spatial variable dictionary and convert the values to arcseconds.
        # =============================================================================
        self.param['horizontal'] = {'spatial': ('r', x)}
        self.param['vertical'] = {'spatial': ('z', y)}

        spatial_name = self.param[self.orientation]['spatial'][0]
        spatial = np.array(self.param[self.orientation]['spatial'][1]) * self.profile.pixscale  # Here is the pixel scale value of the original system
        
        vec = [spatial, sb_r, sb_l, sb, error_up, error_down, c_r, c_l, c, c_err]
        names = [spatial_name, 'sb right', 'sb left', 'sb', 'err up', 'err down', 'ct right', 'ct left', 'ct', 'ct err']
        
        return Table(vec, names=names, meta={'Surface Brighnes': 'table info'})



# =============================================================================
# Shifted radial profiles
# =============================================================================
            
class ProjectedShiftedProfile:
    
    def __init__(self, image, profile, wcs_ori, wcs, zp, pixscale, A=0, rms_sky=0, dx=0, dy=0):
        """
        Class definition to extract a surface brightness profile with an offset from the main axis of the galaxies
        using rectangular apertures over an image. The apertures are built here.
        
        Parameters
        ----------
        image : 2D numpy array
            Data from which to extract the profile.
        Xc : float
            X center coordinate in pixels.
        Yc : float
            Y center coordinate in pixels.
        nbins : int
            Number of bins in the profile.
        npix : int
            Total width of the profile in pixels (number of pixels up to which the profile is extracted).
        height : float
            Height of the rectangular apertures in pixels.
        zp : float
            Zero point [mag/arcsec^2].
        delta : float
            Offset in pixels from the galaxy center at which to extract the profile.
        rms_sky : float
            Sky RMS to include as a source of error.
        A : float
            Galactic absorption or extinction for a given band.
        orientation : str, optional
            Orientation of the profile. 'horizontal' for the galaxy's major axis; 'vertical' for the minor axis.
            Default is 'horizontal'. 

        Returns
        -------
        None.
        """
        
        # Dictionary for all parameters depending on whether the aperture is vertical or horizontal.
        self.param = {'horizontal': {'signo': (-1, 1)},
                      'vertical': {'signo': (1, -1)}}

        # Geometric parameters
        self.orientation = profile.orientation
        
        # WCS
        self.wcs_ori = wcs_ori 
        self.wcs = wcs
        
        # Data parameters
        self.image = image
        self.zp = zp
        self.A = A
        self.rms_sky = rms_sky
        self.pixscale = pixscale
        self.profile = profile

        # =============================================================================
        # Apertures definition
        # =============================================================================
        up_ori = profile.up
        down_ori = profile.down
        
        self.up = ProjectCentralProfile(image, up_ori, wcs_ori, wcs, zp, pixscale, A=A, rms_sky=rms_sky, dx=dx, dy=dy)
        self.down = ProjectCentralProfile(image, down_ori, wcs_ori, wcs, zp, pixscale, A=A, rms_sky=rms_sky, dx=dx, dy=dy)

    def plot(self, color=None, ls='solid'):
        """
        Plot the profile apertures.

        Parameters
        ----------
        color : str, optional
            Colour of the apertures in the plot. Default is None.
        ls : str, optional
            Linestyle ("solid", "dotted", "dashed", "dashdot", etc.). Default is 'solid'.
        cent_ap : bool, optional
            Plot a mark at the central position of the aperture. Default is False.

        Returns
        -------
        Displays a plot of the apertures.
        """
        
        self.up.plot(color=color, ls=ls)
        self.down.plot(color=color, ls=ls)
    
    def saveDS9(self, fname):
        """
        Writes the apertures to a file as DS9 region definitions.

        Parameters
        ----------
        fname : str
            File name.

        Returns
        -------
        A .reg file saved in the current folder.
        """
        self.up.saveDS9(fname + '_up')
        self.down.saveDS9(fname + '_down')        
        

    def surfbrigth(self, method, subpixels=None, sigma=3):
        """
        Extracts the surface brightness profile from the apertures defined in the ShiftedProfile class.

        Parameters
        ----------
        method : str
            Interpolation method as described above.
        subpixels : int, optional
            Number of subpixels if 'subpixel' method is used. Default is None.
        sigma : float, optional
            Sigma value for the sigma-clipping rejection algorithm in the statistics. Default is 3.

        Returns
        -------
        astropy.table.Table
            Table with profile values in counts and surface brightness, and their errors.
        """
        
        # Auxiliary class for converting counts to surface brightness
        c2sb = Counts2Sb(self.zp, self.pixscale, self.A)

        t_1 = self.up.surfbrigth(method, subpixels, sigma)
        t_2 = self.down.surfbrigth(method, subpixels, sigma)
        
        # Retrieve the spatial coordinate name and values directly from one of the tables (they are identical)
        spatial = t_1.columns[0]
        spatial_name = t_1.colnames[0]

        # Right apertures:
        # Counts
        c1_r = t_1['ct right']
        c2_r = t_2['ct right']
        c_r = np.nanmean([c1_r, c2_r], 0)
        
        # Magnitudes for right apertures
        sb_r = c2sb.do(c_r)
        sb_up_r = t_1['sb right']
        sb_down_r = t_2['sb right']
        
        # Left apertures:
        # Counts
        c1_l = t_1['ct left']
        c2_l = t_2['ct left']
        c_l = np.nanmean([c1_l, c2_l], 0)

        # Magnitudes for left apertures
        sb_l = c2sb.do(c_l)
        sb_up_l = t_1['sb left']
        sb_down_l = t_2['sb left']
        
        # Average over all apertures:
        # Counts
        c_mean = np.nanmean([c1_r, c2_r, c1_l, c2_l], 0)
        # Magnitudes for the overall aperture
        sb_mean = c2sb.do(c_mean)
        
        # =============================================================================
        # Errors Calculation
        # =============================================================================
        c1_error = t_1['ct err']
        c2_error = t_2['ct err']        
        
        # Quadratic sum for total error
        err_t = np.sqrt(np.square(c1_error) + np.square(c2_error))
        
        # Error in counts
        c_err = np.copy(err_t)

        up = c_mean + err_t
        down = c_mean - err_t
        
        error_down = -2.5 * np.log10(c_mean / up)
        error_up = -2.5 * np.log10(down / c_mean)

        from astropy.table import Table
        vec = [spatial, sb_r, sb_l, sb_mean, error_up, error_down, c_r, c_l, c_mean, c_err]
        
        names = [spatial_name, 'sb right', 'sb left', 'sb', 'err up', 'err down', 'ct right', 'ct left', 'ct', 'ct err']

        return Table(vec, names=names, meta={'Surface Brighnes': 'table info'})



# ********************************************************************************
# ELLIPTICAL RADIAL PROFILE 
# ********************************************************************************

class EllipticalProfile:
    def __init__(self, data, Xc, Yc, sma, eps, pa, astep=0.1, linear_growth=False, fix_center=False, fix_pa=False, fix_eps=False, ell_threshold=0.1):
        '''
        Class definition to extract elliptical isophote values from an image and obtain a profile.
        This defines the geometry of the initial ellipse.

        Parameters
        ----------
        data : 2D numpy array
            Image containing the data from which to extract the profile.
        Xc : float
            X center coordinate in pixels.
        Yc : float
            Y center coordinate in pixels.
        sma : float
            Semi-major axis of the initial ellipse in pixels (but not the smallest).
        eps : float
            Ellipticity of the initial ellipse.
        pa : float
            Position angle (in radians) of the semi-major axis with respect to the positive x-axis of the image array 
            (rotating towards the positive y-axis). Position angles are defined in the range 0 < PA <= π. 
            Avoid using a starting position angle of 0, since the fitting algorithm may not work properly.
            When ellipses have position angles near either extreme of the range, noise can cause the solution 
            to jump between successive isophotes by nearly 180 degrees.
        astep : float, optional
            Step value for growing/shrinking the semi-major axis. Can be expressed in pixels (if linear_growth=True) 
            or as a relative value (if linear_growth=False). Default is 0.1.
        linear_growth : bool, optional
            Mode for growing/shrinking the semi-major axis. Default is False.
        fix_center : bool, optional
            Fix the ellipse center during fitting? Default is False.
        fix_pa : bool, optional
            Fix the position angle of the ellipse during fitting? Default is False.
        fix_eps : bool, optional
            Fix the ellipticity during fitting? Default is False.
        ell_threshold : float, optional
            Threshold for the object centering algorithm.
            Lowering this value makes the centerer less strict (accepts lower signal-to-noise data). 
            If set very high, the centerer is effectively disabled. In that case, either the supplied geometry is used
            as is, or the fit algorithm will terminate prematurely.
            Default is 0.1.

        Returns
        -------
        None.
        '''
        self.data = data
        self.Xc = Xc
        self.Yc = Yc
        
        # ---------------------------------------------------------------------------
        # Instantiate EllipseGeometry as a global variable of the class
        # ---------------------------------------------------------------------------
        from photutils.isophote import EllipseGeometry     
        self.geometry = EllipseGeometry(self.Xc, self.Yc, sma=sma, eps=eps, pa=pa, astep=astep,
                                        linear_growth=linear_growth, fix_center=fix_center, fix_pa=fix_pa, fix_eps=fix_eps)
        
        # ---------------------------------------------------------------------------
        # Instantiate Ellipse
        # ---------------------------------------------------------------------------
        from photutils.isophote import Ellipse
        self.ellipse = Ellipse(self.data, self.geometry, threshold=ell_threshold)
        
    def find_center(self):
        '''
        Find the center coordinates of the galaxy.

        Returns
        -------
        Upon successful execution, the (x, y) coordinates in the EllipseGeometry instance 
        (i.e., the x0 and y0 attributes) will be modified.
        '''
        self.geometry.find_center(self.data)

    def fit_image(self, sma0=None, minsma=0.0, maxsma=None, step=0.1, conver=0.05, minit=10, maxit=50, fflag=0.5, maxgerr=0.5, sclip=3.0, nclip=3, integrmode='median', linear=None, maxrit=None, fix_center=False, fix_pa=False, fix_eps=False):
        '''
        Fit multiple isophotes to the image array plus an extra fit for the central region.
        Loops over each semi-major axis (sma) length, fitting a single isophote at each sma.
        Returns an IsophoteList instance (a list-like object of Isophote instances sorted by increasing sma).

        Parameters
        ----------
        sma0 : float, optional
            Starting semi-major axis length (pixels). Should be chosen so that the corresponding isophote 
            has good S/N and well-defined geometry. If None, either the sma value from the provided 
            EllipseGeometry is used, or 10 is used as default.
        minsma : float, optional
            Minimum semi-major axis length (pixels). Default is 0.0.
        maxsma : float, optional
            Maximum semi-major axis length (pixels). If None, the algorithm continues until a stopping condition is met.
        step : float, optional
            Step value to grow/shrink the sma (in pixels if linear=True or relative if linear=False). Default is 0.1.
        conver : float, optional
            Convergence criterion. Default is 0.05.
        minit : int, optional
            Minimum number of iterations. Default is 10.
        maxit : int, optional
            Maximum number of iterations. Default is 50.
        fflag : float, optional
            Fitting flag parameter. Default is 0.7.
        maxgerr : float, optional
            Maximum allowed error. Default is 0.5.
        sclip : float, optional
            Sigma clipping value. Default is 3.0.
        nclip : int, optional
            Number of clips. Default is 3.
        integrmode : {'bilinear', 'nearest_neighbor', 'mean', 'median'}, optional
            Integration mode. Default is 'median'.
        linear : bool, optional
            Whether to use linear growth. Default is None.
        maxrit : float, optional
            Maximum iterations. Default is None.
        fix_center : bool, optional
            Whether to fix the center during fitting. Default is False.
        fix_pa : bool, optional
            Whether to fix the position angle during fitting. Default is False.
        fix_eps : bool, optional
            Whether to fix the ellipticity during fitting. Default is False.

        Returns
        -------
        An IsophoteList instance sorted by increasing sma.
        '''
        self.maxsma = maxsma
        self.isolist = self.ellipse.fit_image(sma0, minsma, self.maxsma, step, conver, minit, maxit, fflag, maxgerr, sclip, nclip, integrmode, linear, maxrit, fix_center, fix_pa, fix_eps)
        
        # ---------------------------------------------------------------------------
        # Central ellipse
        # ---------------------------------------------------------------------------
        from photutils.isophote.sample import CentralEllipseSample
        from photutils.isophote.fitter import CentralEllipseFitter
        from photutils.isophote import EllipseGeometry

        # Define central peak position.
        g = EllipseGeometry(self.Xc, self.Yc, 0., 0., 0.)

        # Build a CentralEllipseSample instance and fit the central ellipse.
        sample = CentralEllipseSample(self.data, 0., geometry=g)
        fitter = CentralEllipseFitter(sample)
        center = fitter.fit()
        # Replace the first isophote with the new central fit and sort the list.
        self.isolist[0] = center
        self.isolist.sort()

    def outer_regions(self, sma_cut, step):
        '''
        Perform additional fitting for the outer regions of the object, useful for very large radial distances and/or sky regions.

        Parameters
        ----------
        sma_cut : float
            Radial position (in pixels) from which to start fitting the outer regions.
        step : float
            Step value for growing/shrinking the sma. Use a larger step than in fit_image due to lower S/N.

        Returns
        -------
        An IsophoteList instance that includes the fit from fit_image plus the extended outer regions.
        '''
        from photutils.isophote import EllipseGeometry
        from photutils.isophote import Ellipse
        import numpy as np

        # Find the index of the first isophote with sma greater than sma_cut.
        ind = np.where(self.isolist.sma > sma_cut)[0].min()

        # Get initial parameters from the isophote at the chosen index.
        x0 = self.isolist.x0[ind]
        y0 = self.isolist.y0[ind]
        sma = self.isolist.sma[ind]
        eps = self.isolist.eps[ind]
        pa = self.isolist.pa[ind]

        g = EllipseGeometry(x0, y0, sma, eps, pa)
       
        ellipse = Ellipse(self.data, geometry=g)

        outer_isolist = ellipse.fit_image(integrmode='median', step=step, sma0=sma, minsma=sma, fflag=0.2, sclip=3.0, nclip=3, fix_center=True)

        # Concatenate the outer isophotes with the original list up to index ind.
        self.isolist = self.isolist[:ind] + outer_isolist
        
    def run_over(self, new_data, zp, pixscale, rms_sky=0, A=0, dimming=0, k_corr=0):
        '''
        Applies the previously determined ellipse fitting to a new image using the same isophotes.

        Returns
        -------
        An Astropy Table with surface brightness values and errors.
        '''        
        from photutils.isophote import EllipseSample, Isophote, IsophoteList

        # Loop over all isophotes (excluding the first)
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
        
        # Build an IsophoteList from the new samples.
        new_isolist = IsophoteList(new_isolist)
        
        # ---------------------------------------------------------------------------
        # Central ellipse re-fitting on the new data
        # ---------------------------------------------------------------------------
        from photutils.isophote.sample import CentralEllipseSample
        from photutils.isophote.fitter import CentralEllipseFitter
        from photutils.isophote import EllipseGeometry

        # Define central peak position.
        g = EllipseGeometry(self.Xc, self.Yc, 0., 0., 0.)

        # Build a CentralEllipseSample instance and fit the central ellipse.
        sample = CentralEllipseSample(new_data, 0., geometry=g)
        fitter = CentralEllipseFitter(sample)
        center = fitter.fit()
        new_isolist.append(center)
        new_isolist.sort()
        
        return self.surfbrigth(zp, pixscale, rms_sky=rms_sky, A=A, dimming=dimming, k_corr=k_corr, isolist=new_isolist)

    def surfbrigth(self, zp, pixscale, rms_sky=0, A=0, dimming=0, k_corr=0, isolist=None):
        import numpy as np
        from astropy.table import Table
        # If isolist is not provided, use the original isolist
        if isolist is None:
            isolist = self.isolist
        # Auxiliary class for converting counts to surface brightness        
        c2sb = Counts2Sb(zp, pixscale, A, dimming, k_corr)
        
        # Read counts from the fit
        sma = isolist.sma
        counts = isolist.intens
        counts_err = np.sqrt(np.square(isolist.int_err/np.sqrt(isolist.npix_e)) + np.square(rms_sky/np.sqrt(isolist.npix_e)))
        
        # Signal-to-noise ratio
        SNR = counts / isolist.int_err
        
        # Fitting parameters
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
     
        # Calculate surface brightness with corrections
        sb = c2sb.do(counts)

        # Convert counts errors to magnitude errors
        up = counts + counts_err
        down = counts - counts_err

        error_down = -2.5 * (np.log10(counts / up))
        error_up = -2.5 * (np.log10(down / counts))
        
        # Replace NaN values with 0
        error_up = np.nan_to_num(error_up)
        error_down = np.nan_to_num(error_down)

        # ---------------------------------------------------------------------------
        # Save the results in a table
        # ---------------------------------------------------------------------------
        names = ['sb', 'err up', 'err down', 'ct', 'ct err', 'SNR', 'sma', 'xc', 'xc err', 'yc', 'yc err', 'ell', 'ell err', 'pa', 'pa err', 'b4', 'b4 err']
        vec = [sb, error_up, error_down, counts, counts_err, SNR, sma, xc, xc_err, yc, yc_err, ell, ell_err, pa, pa_err, b4, b4_err]
    
        return Table(vec, names=names, meta={'Surface Brightess': 'table info'})
    
    def plot(self, color=None, ls='solid', smamax=np.inf):
        from photutils import EllipticalAperture
        # ---------------------------------------------------------------------------
        # Read the fitted apertures
        # ---------------------------------------------------------------------------
        print(len(self.isolist))
        x0 = self.isolist.x0
        y0 = self.isolist.y0
        sma = self.isolist.sma
        smna = sma * (1 - self.isolist.eps)
        pa = self.isolist.pa
        # ---------------------------------------------------------------------------
        # Colour palette/map
        # ---------------------------------------------------------------------------
        if color is None:
            from matplotlib.pyplot import cm
            paleta = cm.rainbow(np.linspace(0, 1, len(x0 - 1)))
        else:
            paleta = [color] * (len(x0 - 1))

        for i in range(1, len(x0)):
            if sma[i] <= smamax:
                ap = EllipticalAperture((x0[i], y0[i]), sma[i], smna[i], pa[i])
                ap.plot(color=paleta[i], ls=ls)


# =============================================================================
# Function to correct Halpha images for NII contamination
# =============================================================================

def NIIcorr(flux_Ha, M_B):
    '''
    Function to apply the NII correction to Halpha images.
    See equation B1 in Kennicutt+08.

    Parameters
    ----------
    flux_Ha : np.array
        Array with Halpha + NII flux values already corrected for galactic absorption [erg s^-1 cm^-2].
    M_B : float
        Absolute magnitude in the B band.

    Returns
    -------
    Corrected Halpha flux with NII contamination removed.
    '''
    if M_B <= -21.:
        nIIcorr = 0.54   # log(NII/Halpha) = 0.54
    else:
        nIIcorr = -0.173 * M_B - 3.903
    
    flux_corr = flux_Ha / (1 + 10**nIIcorr)
    
    return flux_corr

# =============================================================================
# Function to convert AB magnitudes to flux in cgs units
# =============================================================================

def magab2fluxcgs(mag_ab, lambda_angs, delta_lambda):
    # Reference: https://en.wikipedia.org/wiki/AB_magnitude 
    flux_jy_nu = 10**(-0.4 * (mag_ab - 8.9))
    flux_cgs_lamda = flux_jy_nu / (3.34e4 * (lambda_angs**2))
    flux_cgs = flux_cgs_lamda * delta_lambda
    
    return flux_cgs

# =============================================================================
# Function to calculate the SFR from Halpha images
# =============================================================================

def SFR_Halpa(flux_cgs, flux_err, DMpc):
    '''
    Function to obtain the Star Formation Rate (SFR) from Halpha data.
    See Kennicutt+2009 and Section 5.7 in Sanchez-Gallego+2012:
        
        SFR(M_sun/yr) = 5.5 × 10^-42 L(Hα) (erg s^-1)

    Parameters
    ----------
    flux_cgs : astropy.table.column.Column
        Column of Halpha flux values in cgs units, corrected for galactic absorption and NII contamination.
    flux_err : astropy.table.column.Column
        Column of the Halpha flux errors in cgs units.
    DMpc : float
        Distance in megaparsecs (Mpc) to the object.

    Returns
    -------
    Astropy Table with SFR values [M_sun/yr] and their errors.
    '''
    
    Flux = flux_cgs
    
    # Convert flux to luminosity: L(Hα) [erg/s]
    # L(Hα) = 4πD^2 * (3.086×10^24)^2 * Flux
    LHa = 4 * np.pi * (DMpc**2) * ((3.086e24)**2) * Flux
    
    # Calculate SFR using the above equation 
    SFR = LHa * 5.5e-42   # Alternatively: 7.9e-42 * LHa (Kennicutt 1998)
        
    # SFR error is proportional to flux error
    SFR_err = (5.5e-42) * 4 * np.pi * (DMpc**2) * ((3.086e24)**2) * flux_err

    vec = [SFR, SFR_err]
    names = ['SFR', 'SFR err']
        
    sfr_table = Table(vec, names=names, meta={'SFR': 'table info'})

    
    # from astropy.table import hstack
    # sfr_table = hstack([SFR, SFR_err])
    # sfr_table['sb'].name = 'SFR'
    # sfr_table['ct err'].name = 'SFR err'
    
    return sfr_table

# =============================================================================
# Function to delete regions from masks interactively
# =============================================================================

def DeleteMaskRegion(mask, image, vmin=0, vmax=0.3):
    """ [Interactive] Function to delete mask regions by clicking on the regions to remove.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from tkinter import messagebox  # tkinter.TkVersion 8.6 
    from astropy.visualization import LogStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    
    plt.ion(), plt.show()

    # Set limits to ensure consistent colour mapping for the masked regions
    vmax_m = np.max(mask)

    remove = True
    plt.figure(1, figsize=(7, 7))
    
    while remove == True:
        happy = False
        while happy == False:
            plt.cla()
            plt.imshow(mask, origin='lower', cmap='tab20c', aspect='equal', vmin=1, vmax=vmax_m, interpolation='none')
            plt.tight_layout()
            plt.show()
        
            mask_new = np.copy(mask)
            image_new = np.copy(image)
            
            # Iteratively take points using ginput
            points = plt.ginput(n=-1, timeout=0, show_clicks=True, mouse_pop=2, mouse_stop=3)
            for point in points:
                coords = np.int_(point)
                value = mask[coords[1]][coords[0]]
                mask_new[mask_new == value] = 0
                image_new[mask_new != 0] = np.nan

            plt.cla()
            # Create normalization object
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch())
            
            plt.imshow(image_new, origin='lower', aspect='equal', norm=norm)
            plt.tight_layout()
            plt.show()
            
            happy = messagebox.askyesno("Masks", "Are you happy with the result?")
            if happy == True:
                mask = np.copy(mask_new)
                remove = messagebox.askyesno("Masks", "Do you want to remove more masked regions?", default='no')
                
    plt.close(1)

    return mask_new


# =============================================================================
# Function to mask a circular section at a given position using the cursor
# =============================================================================

def CircularHandMask(image, mask_value=True, vmin=0, vmax=0.2):
    """ [Interactive] Function to add circular mask regions by specifying the number of new masks and their radius in pixels.
    Then click on the image as many times as required to indicate where to place those masks.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from tkinter import messagebox  # tkinter.TkVersion 8.6 
    from astropy.visualization import LogStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    from photutils import CircularAperture
    plt.ion()
    shape = np.shape(image)
    
    remove = True
    plt.figure(1, figsize=(7, 7))
    
    while remove == True:
        happy = False
        while happy == False:
            plt.cla()
            
            # Create a normalization object
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch())
            plt.imshow(image, origin='lower', aspect='equal', norm=norm)            
            plt.tight_layout()
            plt.show()
        
            image_new = np.copy(image)
            
            # Prompt for mask radius (default is 30 pixels)
            radii = input('masks radii? (default 30 pix): ')
            if radii.isnumeric() == True:
                radii = int(radii)
            else:
                radii = 30
            
            # Interactively get coordinates using ginput
            coords = plt.ginput(n=-1, timeout=0, show_clicks=True, mouse_pop=2, mouse_stop=3)
            
            mask = 0
            for coord in coords:
                # Create a circular aperture at each position
                ap = CircularAperture(coord, radii)
                mask = mask + ap.to_mask(method='center').to_image(shape)
                
            image_new[mask != 0] = mask_value
            
            plt.cla()            
            # Create a normalization object for display
            norm = ImageNormalize(vmin=0., vmax=0.3, stretch=LogStretch())
            
            plt.imshow(image_new, origin='lower', aspect='equal', norm=norm)
            plt.tight_layout()
            plt.show()
            
            happy = messagebox.askyesno("Masks", "Are you happy with the result?")
            if happy == True:
                image = np.copy(image_new)
                remove = messagebox.askyesno("Masks", "Do you want to continue masking this region?")
                
    plt.close(1)

    return image  


def CircularHandMask_new(image, mask_value=True, vmin=0, vmax=0.2):
    """ [Interactive] Function to add circular mask regions by specifying the number of new masks and their radius in pixels.
    Then click on the image as many times as required to indicate where to place those masks.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from tkinter import messagebox  # tkinter.TkVersion 8.6 
    from astropy.visualization import LogStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    from photutils import CircularAperture

    shape = np.shape(image)
    
    remove = True
    plt.figure(1, figsize=(7, 7))
    
    while remove == True:
        happy = False
        while happy == False:
            plt.cla()
            
            # Create a normalization object
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch())
            plt.imshow(image, origin='lower', aspect='equal', norm=norm)            
            plt.tight_layout()
            plt.show()
        
            image_new = np.copy(image)
            im_zeros = np.zeros_like(image)
            
            # Prompt for mask radius (default is 10 pixels)
            radii = input('masks radii? (default 10 pix): ')
            if radii.isnumeric() == True:
                radii = int(radii)
            else:
                radii = 10
            
            # Interactively get coordinates
            coords = plt.ginput(n=-1, timeout=0, show_clicks=True, mouse_pop=2, mouse_stop=3)
            
            mask = 0
            for coord in coords:
                # Create a circular aperture at each position
                ap = CircularAperture(coord, radii)
                mask = mask + ap.to_mask(method='center').to_image(shape)
                
            image_new[mask != 0] = mask_value
            im_zeros[mask != 0] = mask_value
            
            plt.cla()            
            # Create a normalization object for display
            norm = ImageNormalize(vmin=0., vmax=0.3, stretch=LogStretch())
            
            plt.close('all')
            plt.imshow(image_new, origin='lower', aspect='equal', norm=norm)
            plt.tight_layout()
            plt.show()
            
            happy = messagebox.askyesno("Masks", "Are you happy with the result?")
            if happy == True:
                remove = messagebox.askyesno("Masks", "Do you want to continue masking this region?")
                
    plt.close(1)

    return im_zeros  


# =============================================================================
# Elliptical masks
# =============================================================================

def EllipticalHandMask(image, mask_value=True, vmin=0, vmax=0.2):
    """ [Interactive] Function to manually add elliptical masks by specifying how many masks 
    with the same axis values to add.
    Then click on the image as many times as required to indicate where to place those masks.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from tkinter import messagebox  # tkinter.TkVersion 8.6 
    from astropy.visualization import LogStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    from photutils import EllipticalAperture

    remove = True
    plt.figure(1, figsize=(7, 7))
    
    while remove == True:
        happy = False
        while happy == False:
            plt.cla()
            
            # Create a normalization object
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch())
            plt.imshow(image, origin='lower', aspect='equal', norm=norm)            
            plt.tight_layout()
            plt.show()
        
            image_new = np.copy(image)
            
            # Prompt for number of masks to add (default is 1 mask with a=10 and b=5 pixels)
            n = input('How many masks? (default 1): ')
            if n.isnumeric() == True:
                n = int(n)
            else:
                n = 1
                
            a = input('semi-major axis? (default 10): ')
            if a.isnumeric() == True:
                a = int(a)
            else:
                a = 10

            b = input('semi-minor axis? (default 5): ')
            if b.isnumeric() == True:
                b = int(b)
            else:
                b = 5

            # Get coordinates interactively
            coords = plt.ginput(n=n, timeout=0, show_clicks=True)
            
            # Create an elliptical aperture for each selected position
            ap = EllipticalAperture(coords, a=a, b=b)
            
            shape = np.shape(image)
            
            # Loop over each selected position to build the mask
            mask = 0
            for i in range(n):
                mask = mask + ap[i].to_mask(method='center').to_image(shape)
            
            image_new[mask != 0] = mask_value
            
            plt.cla()            
            # Create a normalization object for display
            norm = ImageNormalize(vmin=0., vmax=0.3, stretch=LogStretch())
            
            plt.imshow(image_new, origin='lower', aspect='equal', norm=norm)
            plt.tight_layout()
            plt.show()
            
            happy = messagebox.askyesno("Masks", "Are you happy with the result?")
            if happy == True:
                image = np.copy(image_new)
                remove = messagebox.askyesno("Masks", "Do you want to continue masking this region?")
                
    plt.close(1)

    return image  





# =============================================================================
# Function to partition an image into patches
# =============================================================================

def Patchs(image, size):
    """ Select a patch of an image (data array) of given size in pixels (square region).
    """
    y_l, x_l = np.shape(image)
    # Create sections of the image
    patchs = []
    
    x = np.append(np.arange(0, x_l, size), x_l)
    y = np.append(np.arange(0, y_l, size), y_l)
    
    for i in range(len(x) - 1):
        for j in range(len(y) - 1):
            sec = slice(y[j], y[j+1]), slice(x[i], x[i+1])  
            patchs.append(sec)
    
    return patchs



# =============================================================================
# Truncation fitting
# =============================================================================

def TruncationFitting(x, y, c=1, a1=1, a2=-1, b1=0, sigma=3, c_min=None, c_max=None):
    from astropy.modeling.models import custom_model
    from astropy.modeling import fitting
    from astropy.stats import sigma_clip

    @custom_model
    def TwoLines(x, c=c, a1=a1, a2=a2, b1=b1):
        import numpy as np
        b2 = (a1 - a2) * c + b1
        return np.concatenate([a1 * x[x <= c] + b1, a2 * x[x > c] + b2])

    init = TwoLines()
    init.c.min = c_min
    init.c.max = c_max
    
    fitter = fitting.SLSQPLSQFitter()
    fit = fitting.FittingWithOutlierRemoval(fitter, sigma_clip, niter=3, sigma=sigma)

    return fit(init, x, y)


def TwoLines(x, c=1, a1=1, a2=1, b1=1):
    import numpy as np
    """
    Model of two lines that intersect at c, with slopes a and y-intercept b.
    """
    b2 = (a1 - a2) * c + b1
    return np.concatenate([a1 * x[x <= c] + b1, a2 * x[x > c] + b2])


# =============================================================================
# LogNorm for imshow in log scale
# =============================================================================

def LogNorm(vmin, vmax):    
    from astropy.visualization import LogStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    return ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch())


# =============================================================================
# Functions to convert between kpc and arcsec, and vice versa
# =============================================================================

# Kpc to Arcsec **************************************************************
def Kpc2Arcsec(x_kpc, DMpc):
    """
    Converts a length in kiloparsecs to arcseconds.
    
    Parameters
    ----------
    x_kpc : float
        Value in kiloparsecs.
    DMpc : float
        Distance to the target in megaparsecs.
    """
    import numpy as np
    return np.arctan(x_kpc / (DMpc * 1e3)) * (180. / np.pi) * 60 * 60.


# Arcsec to kpc **************************************************************
def Arcsec2Kpc(x_arcsec, DMpc):
    """
    Converts a length in arcseconds to kiloparsecs.
    
    Parameters
    ----------
    x_arcsec : float
        Value in arcseconds.
    DMpc : float
        Distance to the target in megaparsecs.
    """
    import numpy as np
    return np.tan((x_arcsec / (60. * 60.)) * (np.pi / 180.)) * DMpc * 1e3

          
# Pix to kpc **************************************************************
def Pix2Kpc(pix, scale, DMpc=None, z=None, H0=70, Om0=0.3):
    """
    Converts a length in pixels to kiloparsecs.
    
    Parameters
    ----------
    pix : float
        Value in pixels.
    scale : float
        Pixel scale in arcsec/pixel.
    DMpc : float, optional
        Distance to the target in megaparsecs.
    z : float, optional
        Redshift of the target.
    H0 : float, optional
        Hubble constant. Default is 70.
    Om0 : float, optional
        Matter density parameter. Default is 0.3.
    """
    try:
        import numpy as np
        x_arcsec = pix * scale
        return np.tan((x_arcsec / (60. * 60.)) * (np.pi / 180.)) * DMpc * 1e3
    except:
        # For a given redshift using the angular aperture method:
        from astropy.cosmology import FlatLambdaCDM
        from astropy import units as u
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
        d_A = cosmo.angular_diameter_distance(z=z)
        r_arcsec = (pix * scale) * u.arcsec
        return (r_arcsec * d_A).to(u.kpc, u.dimensionless_angles()).value 


# kpc to pix **************************************************************
def Kpc2Pix(kpc, scale, DMpc=None, z=None, H0=70, Om0=0.3):
    try:
        import numpy as np
        x_arcsec = (np.arctan(kpc / (DMpc * 1e3)) / (np.pi / 180.)) * 60 * 60
        pix = x_arcsec / scale
        return pix
    except:
        from astropy.cosmology import FlatLambdaCDM
        from astropy import units as u
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
        d_A = cosmo.arcsec_per_kpc_proper(z)
        x_arcsec = d_A * kpc * u.kpc 
        return (x_arcsec / scale).value


# =============================================================================
# Standard deviation in a region            
# =============================================================================

class RectStat():
    def __init__(self, data, w, h, n=0, Xc=0, Yc=0, theta=0):
        """ 
        Class to compute statistics over a rectangular region.
        """
        from photutils import RectangularAperture
        
        if n == 0:
            # Geometric parameters of the aperture
            self.ap = [RectangularAperture([Xc, Yc], w=w, h=h, theta=0)]  
        else:
            Xc = np.random.rand(n) * (data.shape[1] - 2 * w) + w
            Yc = np.random.rand(n) * (data.shape[0] - 2 * h) + h
            self.ap = []
            for i in range(n): 
                # Append rectangular apertures with random positions
                self.ap.append(RectangularAperture([Xc[i], Yc[i]], w=w, h=h, theta=0))  
    
        self.counts = np.array([])
        
        for ap in self.ap:
            # Extract the section of data over which to compute statistics
            mask = ap.to_mask(method='center')  # Mask geometry
            mask_data = mask.data  # Mask with zeros and ones based on the method
            sec = mask.cutout(data)  # Section of data matching the mask dimensions
            sec_weight = sec * mask_data  # Weight the data according to the mask
            # Concatenate the data from unmasked pixels
            self.counts = np.concatenate([self.counts, sec_weight[mask_data > 0]])

    def stat(self, sigma_clip=3):
        from astropy.stats import sigma_clipped_stats   
        # Compute mean, median, and standard deviation
        mean, median, stddev = sigma_clipped_stats(self.counts, sigma=sigma_clip, maxiters=10, mask_value=np.nan)
        return [mean, median, stddev]

    def plot(self, color='r', lw=1, ls='solid'):
        for ap in self.ap:
            # Plot each aperture
            ap.plot(color=color, ls=ls, lw=lw)
    
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
# Centering algorithms
# =============================================================================

def Centroid(data, coords, r=3, method='1dg', mask=None):
    """
    Input: image and a list of coordinates. Finds the centroid within a region of radius "r" 
    around the given position in "coords".
    Output: array with the adjusted positions and the distance between the positions "d".

    Parameters
    ----------
    data : numpy array
        Image data.
    coords : numpy array
        Array of coordinates in a (n,2) numpy array.
    r : number, optional
        Size in pixels of the section around the starting position to perform the centroid. Default is 3.
    method : str
        Method for centering from photutils.centroid:
            'com': Center of mass.
            'quadratic': Fit a 2D quadratic to the data.
            '1dg': Fit 1D Gaussians to the marginal distributions.
            '2dg': Fit a 2D Gaussian to the data.
    mask : optional
        Optional mask.

    Returns
    -------
    List containing the adjusted x and y coordinates.
    """    
    
    from photutils.centroids import centroid_com, centroid_quadratic, centroid_1dg, centroid_2dg
    from astropy.stats import sigma_clipped_stats
    import numpy as np
    
    # Dictionary mapping method names to functions
    cent = {'com': centroid_com, 'quadratic': centroid_quadratic, '1dg': centroid_1dg, '2dg': centroid_2dg}
    
    # Initial positions
    px = coords[0]
    py = coords[1]

    # Define the section of the data
    xmas = np.int(px + r)
    xmenos = np.int(px - r)
    ymas = np.int(py + r)
    ymenos = np.int(py - r)
    
    sec = data[ymenos:ymas, xmenos:xmas]

    # Compute the median (sky) from the section and subtract it
    median = sigma_clipped_stats(sec, sigma=3.0)[1]
    sec = sec - median
    
    # Compute the centroid using the specified method
    x_s, y_s = cent[method](sec, mask=mask)
    
    x = x_s + px - r
    y = y_s + py - r
    
    fit_coords = [x, y]
    
    return fit_coords


# =============================================================================
# Photometry along a line
# =============================================================================

class LinePhot:
    def __init__(self, x_a, x_b, y_a, y_b, xmin, xmax, w, h):
        # Line that passes through the two given points
        import numpy as np
        m = (y_b - y_a) / (x_b - x_a)
        y0 = -(y_b - y_a) / (x_b - x_a) * x_a + y_a
        ang = np.arctan(m)
    
        # Define rectangular apertures
        from photutils import RectangularAperture

        delta_x = np.cos(ang)
        delta_y = np.sin(ang)
        
        # Loop to create apertures along the line
        x = xmin
        y = xmin * m + y0
        
        self.ap = []
        
        while x < xmax:
            self.ap.append(RectangularAperture([x, y], w, h, theta=ang))
            x = x + delta_x * w
            y = y + delta_y * w
         
        # Save the parameters for later use
        self.m = m
        self.y0 = y0
        self.ang = ang

    def plot(self, color=None, ls='solid'):
        # Colour palette
        if color is None:
            from matplotlib.pyplot import cm
            paleta = cm.rainbow(np.linspace(0, 1, len(self.ap)))
        else:
            paleta = [color] * (len(self.ap))
        
        for i in self.ap: 
            i.plot(color=paleta[self.ap.index(i)], ls=ls)

    def phot(self, data, method, subpixels=None, sigma=3):
        from astropy.stats import sigma_clipped_stats        
        counts = []
        x = []
        for i in self.ap: 
            sec_data = SecImData(i, data, method, subpixels=subpixels)
            # Compute sigma-clipped median
            median = sigma_clipped_stats(sec_data, sigma=sigma)[1] 
            counts.append(median)
            x.append(i.positions[0])
        return np.array(x), np.array(counts)
    
    def centerLine(self, x):
        return self.m * x + self.y0
    
    def bondary(self, x, px, py):
        m = -1 / self.m
        y0 = py - px * m
        return x * m + y0
    
# =============================================================================
# Sersic + Exponential
# =============================================================================

class SersExp:
    
    def Sersic(self, amplitude, r_eff, n):
        # Instantiate a Sersic1D model with the given parameters.
        from astropy.modeling.models import Sersic1D
        self.Sers = Sersic1D(amplitude=amplitude, r_eff=r_eff, n=n)
        
    def SersFixed(self, amplitude=False, r_eff=False, n=False):
        # Set whether to fix the Sersic parameters during fitting.
        self.Sers.amplitude.fixed = amplitude
        self.Sers.r_eff.fixed = r_eff
        self.Sers.n.fixed = n
        
    def SersBound(self, amplitude=(None, None), r_eff=(None, None), n=(None, None)):
        # Set bounds on the Sersic parameters.
        self.Sers.amplitude.bounds = amplitude
        self.Sers.r_eff.bounds = r_eff
        self.Sers.n.bounds = n

    def Exponential(self, amplitude, tau):
        # Instantiate an Exponential1D model with the given amplitude and scale length (tau).
        from astropy.modeling.models import Exponential1D
        self.Exp = Exponential1D(amplitude=amplitude, tau=tau)
        
    def ExpFixed(self, amplitude=False, tau=False):
        # Set whether to fix the Exponential parameters during fitting.
        self.Exp.amplitude.fixed = amplitude
        self.Exp.tau.fixed = tau
        
    def ExpBound(self, amplitude=(None, None), tau=(None, None)):
        # Set bounds on the Exponential parameters.
        self.Exp.amplitude.bounds = amplitude
        self.Exp.tau.bounds = tau

    # The following Fitter method is commented out.
    # def Fitter(self, fitter=None):
    #     if fitter is None:
    #         from astropy.modeling.fitting import SLSQPLSQFitter 
    #         self.fitter = SLSQPLSQFitter()
    #     else:
    #         from astropy.modeling.fitting import fitter 
    #         self.fitter = fitter()

    def Fit(self, r, ct, weights=None, sigma=3):
        # Fit the combined model (Sersic + Exponential) to the data.
        from astropy.stats import sigma_clip
        from astropy.modeling.fitting import FittingWithOutlierRemoval
        from astropy.modeling.fitting import SLSQPLSQFitter 

        init = self.Sers + self.Exp  
        fitter = SLSQPLSQFitter()
        fit = FittingWithOutlierRemoval(fitter, sigma_clip, niter=3, sigma=sigma)
        return fit(init, r, ct, weights=weights)

    
# =============================================================================
# HI Surface Brightness (Surface Mass Density) Profiles from VLA Data
# =============================================================================

def surfDensVLAHI(data, bmaj, bmin, freq):
    '''
    Function to compute the surface mass density of HI from VLA intensity maps in Jy/beam km/s.
    
    First, converts between surface brightness (flux density per area) and brightness temperature.
    Conversion from Jy/beam to K is performed by specifying the beam area.
    See:
      https://docs.astropy.org/en/stable/units/equivalencies.html
      https://docs.astropy.org/en/stable/api/astropy.units.equivalencies.brightness_temperature.html

    Then converts (K km/s) to Msun/pc^2.

    Parameters
    ----------
    data : astropy.table.column.Column
        HI intensity values from VLA maps in Jy/beam km/s.
    bmaj : float
        Beam major axis (BMAJ) in degrees.
    bmin : float
        Beam minor axis (BMIN) in degrees.
    freq : float
        Rest frequency (RESTFREQ) in Hz.

    Returns
    -------
    Surface mass density in Msun/pc^2 as an astropy.table.column.Column.
    '''
    
    import numpy as np
    from astropy import units as u
    
    bmaj = bmaj.to(u.arcsec)
    bmin = bmin.to(u.arcsec)
    fwhm_to_sigma = 1. / (8 * np.log(2))**0.5
    beam_area = 2. * np.pi * (bmaj * bmin * fwhm_to_sigma**2)
    '''
    The following alternative method is commented out:
    
    freq = freq.to(u.GHz)
    equiv = u.brightness_temperature(freq)
    brigTemp = (data * u.Jy / beam_area).to(u.K, equivalencies=equiv)
    n_HI = 1.823e18 * brigTemp / u.K  # HI column density [atoms/cm^2]
    atomHI = 1.673534e-27 * u.kg       # Mass of a hydrogen atom
    surfMassDens = (n_HI * atomHI / (u.cm * u.cm))
    surfMassDens = surfMassDens.to(u.Msun / (u.parsec * u.parsec))
    '''
    
    # Alternative formula (eq. 2 in https://arxiv.org/pdf/2110.01618.pdf)
    # Gives the same result as above WITHOUT inclination correction.
    surfMassDens = 8794 * data / (bmaj * bmin)    
    
    # If correcting for inclination, one could use:
    # surfMassDens = 8794 * data * np.cos(88.5 * np.pi / 180.) / (bmaj * bmin)

    return surfMassDens     


# =============================================================================
# Look for GAMA groups within an RA, DEC region
# =============================================================================

def table_gama_reg(ramin, ramax, decmin, decmax, tablepath):
    '''
    Retrieves GAMA groups within a specified region of the sky.

    Parameters
    ----------
    ramin : float
        Minimum RA value of the region in degrees.
    ramax : float
        Maximum RA value of the region in degrees.
    decmin : float
        Minimum DEC value of the region in degrees.
    decmax : float
        Maximum DEC value of the region in degrees.
    tablepath : str
        Path to the GAMA groups table.

    Returns
    -------
    None.
    '''
    # (Function body not implemented.)
    pass


# =============================================================================
# Circular Profiles
# =============================================================================

class CircularApProfile:
    """ 
    Class to create circular aperture profiles.
    """
    def __init__(self, image, Xc, Yc, nbins, r_min, npix):        
        self.image = image
        self.Xc = Xc
        self.Yc = Yc
        self.nbins = nbins
        self.npix = npix
        self.r_min = r_min 
        
        # Define the parameters for the apertures.
        import numpy as np
        i = np.log10(r_min)
        l = np.log10(self.npix)
        step = np.logspace(i, l, self.nbins)
        self.r = step
         
        # Build apertures: the first aperture is a circle and the rest are concentric annuli.
        from photutils import CircularAnnulus
        self.ap = []
        for i in range(len(self.r) - 1):
            self.ap.append(CircularAnnulus((self.Xc, self.Yc), r_in=self.r[i], r_out=self.r[i+1]))

    def plot(self, image=None, color=None, ls='solid', cmap='viridis', vmin=None, vmax=None):
        # Colour palette.
        if color is None:
            from matplotlib.pyplot import cm
            paleta = cm.rainbow(np.linspace(0, 1, len(self.r)))
        else:
            paleta = [color] * (len(self.ap))
         for i in range(len(self.ap)):
             self.ap[i].plot(color=paleta[i], ls=ls)
     
     def sbprofile(self, zp, pixscale, method='center', rms_sky=0, A=0, dimming=0, k_corr=0, subpixels=None, sigma=3):
         from astropy.stats import sigma_clipped_stats
         # First, calculate the median within each annulus.
         cts = []
         c_err = []
         for i in self.ap:
             sec_data = SecImData(i, self.image, method, subpixels=subpixels)
             median = sigma_clipped_stats(sec_data, sigma=sigma)[1]
             cts.append(median)

             # Errors:
             # Poisson error
             error_poiss = np.sqrt(np.nanmean(np.square(median - sec_data)) / np.size(sec_data))
             # Bin error from sky RMS
             error_bin = rms_sky / np.sqrt(float(np.size(sec_data)))
             # Total error
             err_t = np.sqrt(np.square(error_poiss) + np.square(error_bin))
             c_err.append(err_t)
         
         c_err = np.array(c_err)
         cts = np.array(cts)
 
         # Compute magnitude errors
         up = cts + c_err
         down = cts - c_err
         error_down = -2.5 * np.log10(cts / up)
         error_up = -2.5 * np.log10(down / cts)

         # Define "r" as the midpoint of each annulus.
         c2sb = Counts2Sb(zp, pixscale, A=0, dimming=0, k_corr=0)
         r = self.r[:-1] + np.diff(self.r) / 2
         sb = c2sb.do(cts)
         
         vec = [r, sb, error_up, error_down, cts, c_err]
         names = ['r', 'sb', 'err up', 'err down', 'ct', 'ct err']
         from astropy.table import Table
         return Table(vec, names=names, meta={'Surface Brightnes': 'table info'})
     
     def model2d(self, method, subpixels=None, sigma=3):
         from astropy.stats import sigma_clipped_stats
         # First, calculate the median within each annulus.
         cts = []
         for i in self.ap:
             sec_data = SecImData(i, self.image, method, subpixels=subpixels)
             median = sigma_clipped_stats(sec_data, sigma=sigma)[1]
             cts.append(median)

         # Define "r" as the midpoint of each annulus.
         r = np.insert(self.r[1:] + np.diff(self.r) / 2, 0, self.r[0] / 2)
             
         # Build the 2D model by interpolating the intensity profile.
         from scipy.interpolate import LSQUnivariateSpline
         from photutils.isophote import EllipseGeometry
         finely_spaced_sma = np.arange(r[0], r[-1], 0.1)
         nodes = r[2:-2]
         shape = self.image.shape
         intens_array = LSQUnivariateSpline(r, cts, nodes)(finely_spaced_sma)
         x0 = self.Xc
         y0 = self.Yc

         result = np.zeros(shape=shape)
         weight = np.zeros(shape=shape)

         # For each interpolated ring, assign intensity values to the output image.
         for index in range(1, len(finely_spaced_sma)):
             sma0 = finely_spaced_sma[index]
             geometry = EllipseGeometry(x0, y0, sma0, eps=0, pa=0)
             intens = intens_array[index]
    
             # Scan angles slightly beyond 360 degrees to ensure full coverage.
             r_val = sma0
             phi = 0.
             while phi <= 2 * np.pi + 0.05:
                 harm = 0.
                 x = r_val * np.cos(phi) + x0
                 y = r_val * np.sin(phi) + y0
                 i_pix = int(x)
                 j_pix = int(y)
    
                 if (i_pix > 0 and i_pix < shape[1] - 1 and j_pix > 0 and j_pix < shape[0] - 1):
                     fx = x - float(i_pix)
                     fy = y - float(j_pix)
    
                     result[j_pix, i_pix] += (intens + harm) * (1. - fy) * (1. - fx)
                     result[j_pix, i_pix + 1] += (intens + harm) * (1. - fy) * fx
                     result[j_pix + 1, i_pix] += (intens + harm) * fy * (1. - fx)
                     result[j_pix + 1, i_pix + 1] += (intens + harm) * fy * fx
    
                     weight[j_pix, i_pix] += (1. - fy) * (1. - fx)
                     weight[j_pix, i_pix + 1] += (1. - fy) * fx
                     weight[j_pix + 1, i_pix] += fy * (1. - fx)
                     weight[j_pix + 1, i_pix + 1] += fy * fx
    
                     phi = max((phi + 0.75 / r_val), 0.05)
                     r_val = max(geometry.radius(phi), 0.5)
                 else:
                     phi = max((phi + 0.75 / r_val), 0.05)
                     r_val = max(geometry.radius(phi), 0.5)
    
         weight[np.where(weight <= 0.)] = 1.
         result /= weight
         return result


# =============================================================================
# ELLIPTICAL APERTURES RADIAL PROFILE 
# =============================================================================
class EllipticalApProfile:
    """ 
    Class to extract an elliptical radial profile using elliptical apertures.
    """
    def __init__(self, image, Xc, Yc, sma_arr, ell_arr, pa_arr):
        self.image = image
        self.Xc = Xc
        self.Yc = Yc
        self.ell = ell_arr
        self.theta = pa_arr
        self.a = sma_arr
        
        # Define the semi-minor axes for the ellipses.
        self.b = (1 - self.ell) * self.a

        # Build apertures: the first aperture is an ellipse and the rest are concentric annuli.
        from photutils import EllipticalAperture, EllipticalAnnulus
        self.ap = [EllipticalAperture((self.Xc, self.Yc), a=self.a[0], b=self.b[0], theta=self.theta[0])]
        for i in range(len(self.a) - 1):
            self.ap.append(EllipticalAnnulus((self.Xc, self.Yc), a_in=self.a[i], a_out=self.a[i+1], b_out=self.b[i+1], theta=self.theta[i]))

    def plot(self, image=None, color=None, ls='solid', cmap='viridis', vmin=None, vmax=None):
        import matplotlib.pyplot as plt
        # If no image is provided, use the original image.
        if np.shape(image) == ():
            image = self.image

        # Display the image.
        plt.imshow(image, cmap=cmap, vmin=vmin, vmax=vmax, origin='lower')
        
        # Colour palette for apertures.
        if color is None:
            from matplotlib.pyplot import cm
            paleta = cm.rainbow(np.linspace(0, 1, len(self.a)))
        else:
            paleta = [color] * (len(self.a))
            
        for i in range(len(self.ap)):
            self.ap[i].plot(color=paleta[i], ls=ls)
        
    def surfbrigth(self, zp, pixscale, method='center', rms_sky=0, A=0, dimming=0, k_corr=0, subpixels=None, sigma=3):
        from astropy.table import Table
        from astropy.stats import sigma_clipped_stats
        
        self.zp = zp
        self.rms_sky = rms_sky
        self.A = A
        self.pixelscale = pixscale
        self.dimming = dimming
        self.k_corr = k_corr
        
        # The first element is the central ellipse; the rest are annuli.
        cts = []
        c_err = []
        sb = []
        error_up = []
        error_down = []
        noise = []

        for i in range(len(self.ap)):
            mask = self.ap[i].to_mask(method=method, subpixels=subpixels)
            mask_data = mask.data
            sec = mask.cutout(self.image)
            sec_weight = sec * mask_data
            sec_data = sec_weight[mask_data > 0]
            
            # Intensity in the bin.
            median = sigma_clipped_stats(sec_data, sigma=sigma)[1]
            cts.append(median)

            # Errors:
            # Poisson error.
            error_poiss = np.sqrt(np.nanmean(np.square(median - sec_data)) / float(np.size(sec_data)))
            # Bin error from sky RMS.
            error_bin = self.rms_sky / np.sqrt(float(np.size(sec_data)))
            noise.append(error_bin)
        
            # Total error in counts.
            err_t = np.sqrt(np.square(error_poiss) + np.square(error_bin))
            c_err.append(err_t)
            
        c_err = np.array(c_err)
        cts = np.array(cts)

        # Magnitude errors.
        up = cts + c_err
        down = cts - c_err
        error_down = -2.5 * np.log10(cts / up)
        error_up = -2.5 * np.log10(down / cts)

        # Define "r" as the midpoint of each annulus.
        r = np.insert(self.a[1:] + np.diff(self.a) / 2, 0, self.a[0] / 2)
        c2sb = Counts2Sb(self.zp, self.pixelscale, A=self.A, dimming=self.dimming, k_corr=self.k_corr)
        sb = c2sb.do(cts)
        
        # Signal-to-noise ratio.
        SNR = cts / noise
        
        vec = [r, sb, error_up, error_down, cts, c_err, SNR]
        names = ['r', 'sb', 'err up', 'err down', 'ct', 'ct err', 'SNR']
           
        return Table(vec, names=names, meta={'Surface Brightnes': 'table info'})
    
# =============================================================================
# Scale factor for subtracting a PSF scattered field of stars
# =============================================================================

class ScalePSF:
    """ 
    Class to calculate a scaling factor for subtracting the scattered PSF light from stars.
    """
    def __init__(self, image, psf, star_center, f_min, f_max, r_min=50, r_max=400, psf_mask=None, psf_center='mid', step=1):
        self.image = image
        self.psf = psf
        self.f_min = f_min
        self.f_max = f_max
        self.star_center = star_center
        
        # Determine the center of the PSF.
        if psf_center == 'mid':
            psf_center = [int(psf.shape[0] / 2) + 1, int(psf.shape[1] / 2) + 1]
        
        # Define the parameters for the apertures.
        i = np.log10(r_min)
        l = np.log10(r_max + step)
        bins = int((r_max + step - r_min) / step)
        r = np.logspace(i, l, bins)

        # Build apertures for the star image and PSF image.
        from photutils import CircularAnnulus
        self.img_ap = []
        self.psf_ap = []
        for i in range(len(r) - 1):
            self.img_ap.append(CircularAnnulus(star_center, r_in=r[i], r_out=r[i+1]))
            self.psf_ap.append(CircularAnnulus(psf_center, r_in=r[i], r_out=r[i+1]))

        # Extract profiles from the apertures.
        img_profile = []
        psf_profile = []
        psf_masked = np.copy(psf)
        psf_masked[psf_mask != 0] = np.nan
            
        for i, j in zip(self.img_ap, self.psf_ap):
            # For the star image.
            sec_data = SecImData(i, self.image, method='center')
            median = sigma_clipped_stats(sec_data, sigma=3)[1]
            img_profile.append(median)

            # For the PSF image.
            sec_data = SecImData(j, psf_masked, method='center')
            median = sigma_clipped_stats(sec_data, sigma=3)[1]
            psf_profile.append(median)
        
        self.r = r[:-1] + np.diff(r) / 2
        self.img_profile = np.array(img_profile)
        self.psf_profile = np.array(psf_profile)
        self.selected = (self.img_profile < self.f_max) & (self.img_profile > self.f_min)

    def model_fluxRange(self, model=None):
        # Calculate scaling factor based on the ratio of the star and PSF profiles within a specified range.
        self.scale_factor_dist = np.log10(self.img_profile[self.selected] / self.psf_profile[self.selected])
        self.scale_factor = sigma_clipped_stats(self.scale_factor_dist, cenfunc='median', sigma=2.5)[0]
    
        if model is None:
            model = np.zeros(self.image.shape)
            
        for i_psf in np.arange(self.psf.shape[0]):
            for j_psf in np.arange(self.psf.shape[1]):
                i = i_psf + int(self.star_center[1]) - int(self.psf.shape[1] / 2)
                j = j_psf + int(self.star_center[0]) - int(self.psf.shape[0] / 2)
                if (i > 0 and i < self.image.shape[1] - 1 and j > 0 and j < self.image.shape[0] - 1 and not np.isnan(self.psf[i_psf, j_psf])):
                    model[i, j] = self.psf[i_psf, j_psf] * 10**self.scale_factor
        return model

    def model_optimize(self, scale_factor=0, n_bs=100, model=None):
        # Optimization of the scaling factor.
        self.scale_factor_dist = []
        if scale_factor == 0:
            for _ in range(n_bs):
                l = len(self.img_profile[self.selected])
                s = int(l * 0.68)
                index = np.random.randint(0, l, size=s)
                def objective_function(a):
                    img_profile = self.img_profile[self.selected][index]
                    psf_profile = self.psf_profile[self.selected][index]
                    return np.nansum(np.square((img_profile - a * psf_profile) * self.r[self.selected][index]))
                res = minimize_scalar(objective_function)
                self.scale_factor_dist.append(res.x)
            self.scale_factor = np.median(self.scale_factor_dist)
            print('Optimized sf :' + str(self.scale_factor))
        if scale_factor != 0:
            self.scale_factor = scale_factor
            print('Manual sf :' + str(self.scale_factor))
        if model is None:
            model = np.zeros(self.image.shape)
        # Build the 2D model using the optimized scaling factor.
        for i_psf in np.arange(self.psf.shape[0]):
            for j_psf in np.arange(self.psf.shape[1]):
                i = i_psf + int(self.star_center[1]) - int(self.psf.shape[1] / 2)
                j = j_psf + int(self.star_center[0]) - int(self.psf.shape[0] / 2)
                if (i > 0 and i < self.image.shape[1] - 1 and j > 0 and j < self.image.shape[0] - 1 and not np.isnan(self.psf[i_psf, j_psf])):
                    model[i, j] = self.psf[i_psf, j_psf] * self.scale_factor
        return model

    def plot(self, savefig=None):
        fig, ax = plt.subplots(2, 2, figsize=(6, 6))
        
        # Star profile:
        coutout = self.img_ap[-1].to_mask(method='exact').cutout(self.image)
        ax[0][0].imshow(coutout, origin='lower', cmap='viridis', norm=LogNorm(vmin=0, vmax=10))
        x_c, y_c = coutout.shape[0] / 2, coutout.shape[1] / 2

        # Plot inner and outer annuli for the star profile.
        circ_int = CircularAperture([x_c, y_c], self.img_ap[0].r_in)
        circ_int.plot(axes=ax[0][0], ls='--', color='r')
    
        circ_int = CircularAperture([x_c, y_c], self.img_ap[-1].r_out)
        circ_int.plot(axes=ax[0][0], ls='--', color='r')
        ax[0][0].set_xticks([])
        ax[0][0].set_yticks([])
        
        # PSF profile:
        coutout = self.psf_ap[-1].to_mask(method='exact').cutout(self.psf)
        ax[0][1].imshow(coutout, origin='lower', cmap='viridis', norm=LogNorm(vmin=0, vmax=0.00001))
        x_c, y_c = coutout.shape[0] / 2, coutout.shape[1] / 2
        circ_int = CircularAperture([x_c, y_c], self.psf_ap[0].r_in)
        circ_int.plot(axes=ax[0][1], ls='--', color='r')
    
        circ_int = CircularAperture([x_c, y_c], self.psf_ap[-1].r_out)
        circ_int.plot(axes=ax[0][1], ls='--', color='r')
        ax[0][1].set_xticks([])
        ax[0][1].set_yticks([])
        
        # Plot the profiles:
        ax[1][0].plot(self.r, self.img_profile, '-o', markersize=2, label='Star profile')
        ax[1][0].plot(self.r, self.psf_profile * self.scale_factor, '--o', markersize=2, label='PSF profile')
        ax[1][0].axhline(self.f_max, ls='--', color='k')
        ax[1][0].axhline(self.f_min, ls='--', color='k')
        ax[1][0].legend()
        ax[1][0].set_xscale('log')
        ax[1][0].set_yscale('log')
        ax[1][0].set_xlim([1, self.img_ap[-1].r_out + 30])
       
        # Histogram of the scaling factor distribution:
        ax[1][1].hist(self.scale_factor_dist, bins='auto')
        ax[1][1].axvline(self.scale_factor, color='r', ls='--')
        
        plt.tight_layout()
        plt.show()  

        if savefig is not None:
            plt.savefig(savefig)


# =============================================================================
# Function: returns distances map from the centre of each galaxy [Mireia]
# =============================================================================

def dist_ellipse_map(N, xc, yc, ratio, pos_ang):  # MIREIA's!!
    import numpy as np
    # Convert pos_ang from degrees to radians.
    ang = pos_ang * np.pi / 180.
    cosang = np.cos(ang)
    sinang = np.sin(ang)
    if len(N) == 2:
        nx = N[0]
        ny = N[1]
    elif len(N) == 1:
        ny = N[0]
        nx = N[0]
    # Create x and y coordinate arrays, shifting so that the galaxy centre is at (xc, yc)
    x = np.arange(0, nx, 1) - [xc - 1.] * nx
    y = np.arange(0, ny, 1) - [yc - 1.] * ny
    im = np.zeros([ny, nx])
    # Rotate the coordinate grid to align with the ellipse orientation.
    xcosang = x * cosang
    xsinang = x * sinang
    for ii in np.arange(0, ny, 1):
         xtemp = xcosang + y[ii] * sinang
         ytemp = -xsinang + y[ii] * cosang
         # The ratio corresponds to 1 - ellipticity, i.e., b/a.
         # Equation of an ellipse: 1 = (x/a)^2 + (y/b)^2
         im[ii, :] = np.sqrt((ytemp / ratio) ** 2 + xtemp ** 2)  # Semi-major axis (SMA)
    return im


def dist_galaxies(image, x_gals, y_gals, ratio_gals, pa_gals):
    # Get the dimensions of the image.
    ny, nx = np.shape(image)
    dist = []
    # For each galaxy, compute its distance map using its parameters.
    for tt in np.arange(0, len(x_gals)):
        grid = dist_ellipse_map([nx, ny], x_gals[tt], y_gals[tt], ratio_gals[tt], pa_gals[tt])
        grid = np.array(grid, dtype=int)
        dist.append(grid)
    # Return the pixel-wise minimum distance from all galaxies.
    r = np.amin(dist, axis=0)
    return r  # 2D array of minimum distances (r) to the centre of galaxies


# =============================================================================
# Function: extracts the (counts) profile using the distances map
# =============================================================================

def dist_ellipse_prof(masked_data, dist_map, sma0, step, nbins, SkyRMS, zp, pixscale, A=0, dimming=0, k_corr=0):
    '''
    Extracts a radial profile from masked image data using a distances map.

    Parameters
    ----------
    masked_data : numpy array 
        Image values from the masked_data FITS file.
    dist_map : numpy array 
        Distances map (e.g., from the galaxy centre) from the dist_map FITS file.
    sma0 : int
        Starting value for the first aperture length in pixels.
    step : float
        Step value for increasing the aperture size (in pixels).
    nbins : int
        Number of apertures for the profile.
    SkyRMS : float
        Sky background RMS.
    zp : float
        Zero point [mag].
    pixscale : float
        Pixel scale in arcsec/pixel.
    A : float, optional
        Galactic absorption/extinction (default 0).
    dimming : float, optional
        Dimming correction in mag/arcsec^2 (default 0).
    k_corr : float, optional
        K-correction (default 0).

    Returns
    -------
    astropy.table.Table
        Table with surface brightness values and errors.
    '''
    from astropy.table import Table
    from astropy.stats import sigma_clipped_stats

    ell_bin = [sma0 * (1. + step)]
    ct = []
    r = []
    c_err = []
    noise = []
    
    for nbin in range(nbins):
        # Compute the boundaries of the current bin.
        ell_bin.append(ell_bin[nbin] * (1. + step))
        # Extract data within the current annulus defined by the distance map.
        data_bin = masked_data[(dist_map >= ell_bin[nbin]) & (dist_map < ell_bin[nbin + 1])]
        
        # Compute the sigma-clipped median of the bin.
        median = sigma_clipped_stats(data_bin, sigma=1, mask_value=np.nan)[1]
        ct.append(median)
        
        # Compute the mid-point of the current annulus.
        r.append(ell_bin[nbin] + (ell_bin[nbin + 1] - ell_bin[nbin]) / 2)
    
        # =============================================================================
        # Errors:
        # =============================================================================
        # Poisson error.
        error_poiss = np.sqrt(np.nanmean(np.square(median - data_bin)) / float(np.size(data_bin)))
        # Bin error from sky RMS.
        error_bin = SkyRMS / np.sqrt(float(np.size(data_bin)))
        noise.append(error_bin)
    
        # Total error.
        err_t = np.sqrt(np.square(error_poiss) + np.square(error_bin))
        c_err.append(err_t)

    # Convert lists to arrays.
    r = np.array(r)
    ct = np.array(ct)
    ct_err = np.array(c_err)

    # Calculate signal-to-noise ratio.
    SNR = ct / noise
 
    # Convert counts to surface brightness.
    c2sb = Counts2Sb(zp, pixscale, A=A, dimming=dimming, k_corr=k_corr)
    sb = c2sb.do(ct)
    sb_err = c2sb.do_errors(ct_err)
    
    # Build table.
    vec = [r, ct, ct_err, SNR, sb, sb_err[0], sb_err[1]]
    names = ['r', 'ct', 'ct err', 'SNR', 'sb', 'err up', 'err down']
       
    return Table(vec, names=names, meta={'Counts and SB Profile': 'table info'})
 
 
# =============================================================================
# Function: extracts the radius at which the ICL dominates the flux from the distances map
# =============================================================================

def findRadius(image, center=None):
     # For center, provide [x, y]
     imsize = np.shape(image)

     if center is None:
         center = [int(imsize[1] / 2.), int(imsize[0] / 2.)]

     distIm = dist_ellipse_map(imsize, center[0], center[1], 1, 0)
     d_vec = np.arange(0, 2000, 10)
     profiles = []
     radius = []

     for nn in np.arange(0, len(d_vec)):
         ind = (distIm < d_vec[nn + 1]) & (distIm > d_vec[nn])
         median_value = np.nanmedian(image[ind])

         if (len(image[ind][np.isnan(image[ind])]) > (len(image[ind]) / 2)):
             median_value = np.nan
         profiles.append(median_value)
         radius.append(d_vec[nn] + (d_vec[nn + 1] - d_vec[nn]) / 2.)
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
    
    # For each isophote value:
    for i, isoph_value in enumerate(isoph_values):
        # 1. Calculate median counts value.
        isoph_pixs_val.append(masked_data[isophotes_im == isoph_value])
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
    
        # Errors:
        error_poiss = np.sqrt(np.nanmean(np.square(ct_median - isoph_pixs_val[i])) / float(np.size(isoph_pixs_val[i])))
        error_bin = SkyRMS / np.sqrt(float(np.size(isoph_pixs_val[i])))
        noise.append(error_bin)
    
        err_t = np.sqrt(np.square(error_poiss) + np.square(error_bin))
        c_err.append(err_t)

    # Convert lists to arrays.
    r = np.array(r)
    ct = np.array(ct)
    ct_err = np.array(c_err)

    # Compute signal-to-noise ratio.
    SNR = ct / noise
 
    # Convert counts to surface brightness.
    c2sb = Counts2Sb(zp, pixscale, A=A, dimming=dimming, k_corr=k_corr)
    sb = c2sb.do(ct)
    sb_err = c2sb.do_errors(ct_err)
    
    # Build table.
    vec = [r, ct, ct_err, SNR, sb, sb_err[0], sb_err[1]]
    names = ['r', 'ct', 'ct err', 'SNR', 'sb', 'err up', 'err down']
       
    return Table(vec, names=names, meta={'Counts and SB Profile': 'table info'})
 
 
# =============================================================================
# Function to save images in .fits format   
# =============================================================================
 
def save_fits_image(data, header, save_path, new_file_name):
    """
    Function to write and save a FITS image.

    Parameters
    ----------
    data : numpy array
        Image data.
    header : FITS header
        Header information.
    save_path : str
        Path to save the new file.
    new_file_name : str
        Name of the new FITS file.

    Returns
    -------
    None.
    """
    from astropy.io import fits
    from astropy.wcs import WCS
    # Create a PrimaryHDU object to encapsulate the data.
    hdu = fits.PrimaryHDU(data)
    wcs = WCS(header)
    hdu.header.update(wcs.to_header())
    hdu.writeto(save_path + new_file_name, overwrite=True)  # Write the file, overwriting if necessary.


def flux_in_ellipse(data, Xc, Yc, maxsma, ell, pa):
    '''
    Estimates the total flux (in counts) within an elliptical aperture.

    Parameters
    ----------
    data : numpy array
        Image data.
    Xc : float
        X coordinate (in pixels) of the ellipse centre.
    Yc : float
        Y coordinate (in pixels) of the ellipse centre.
    maxsma : float
        Semi-major axis length of the ellipse (in pixels).
    ell : float
        Ellipticity of the ellipse.
    pa : float
        Position angle of the ellipse.
    rms_sky : float, optional
        Sky RMS (not used in current implementation).

    Returns
    -------
    Tuple containing:
        - Total flux within the ellipse (in counts)
        - Associated error
        - The elliptical aperture object
    '''
    from photutils import EllipticalAperture
    b = (1 - ell) * maxsma
    Ell_Ap = EllipticalAperture((Xc, Yc), a=maxsma, b=b, theta=ell)
    Ell_Section = SecImData(Ell_Ap, data, method='center', subpixels=None)
        
    flux = np.nansum(Ell_Section)
    # Compute error as the standard deviation of the values within the aperture.
    error = np.sqrt(np.std(Ell_Section))
    return flux, error, Ell_Ap


def patch_mask_cold(im_chunk):
    '''
    Function to generate a mask for a given image chunk as part of a multiprocessing algorithm.
    (Not finished; needs adjustments.)

    Parameters
    ----------
    im_chunk : numpy array
        A chunk of the image.

    Returns
    -------
    mask_chunk : numpy array
        Mask for the image chunk.
    '''
    # The following code is commented out and not implemented:
    # threshold = detect_threshold(im_chunk, nsigma=1.1, sigclip_sigma=2, mask_value=np.nan)
    # mask_chunk = np.zeros_like(im_chunk)
    # try:
    #     segm = detect_sources(im_chunk, threshold, npixels=1200, filter_kernel=kernel)
    #     mask_chunk = deblend_sources(im_chunk, segm, npixels=1200, filter_kernel=kernel).data
    # except:
    #     print('No object detected in patch', patch)
    # return mask_chunk


def magnitudes(ct_array, zp, DMpc=0, z=0, k_corr=0, band='i'):
    '''
    Calculates the absolute magnitude of a galaxy/structure from flux values in counts.

    Steps:
        1. Compute apparent magnitude: m = -2.5*log10(ct_array) + zp - k_corr.
        2. Compute absolute magnitude: M = m - 5*log10(distance in pc).

    Solar absolute magnitudes (AB system): {'g': 5.11, 'r': 4.65, 'i': 4.53}
    (Reference: Table 3 in Willmer2018)

    Parameters
    ----------
    ct_array : numpy array
        Flux values in counts.
    zp : float
        Zero point in magnitudes.
    DMpc : float, optional
        Distance to the object in Mpc.
    z : float, optional
        Redshift of the source.
    k_corr : float, optional
        K-correction.
    band : str, optional
        Photometric band (default 'i').

    Returns
    -------
    dict
        Dictionary with keys for apparent ('m_band') and absolute ('M_band') magnitudes.
    '''
    import numpy as np
    mag_array = -2.5 * np.log10(ct_array) + zp - k_corr

    # For nearby sources:
    if DMpc != 0:    
        abs_mag = mag_array - 5 * np.log10((DMpc * 1e6) / 10)
    
    # For sources at a redshift (using luminosity distance)
    if z != 0:
        from astropy.cosmology import FlatLambdaCDM
        import astropy.units as u
        cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)
        dist_lum = cosmo.luminosity_distance(0.2)
        abs_mag = mag_array - 5 * np.log10((dist_lum.value * 1e6) / 10)
        
    m = 'm_' + band
    M = 'M_' + band
    
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
    abs_mag : float
        Absolute magnitude of the object.
    band : str, optional
        Photometric band (default 'i').

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
        Luminosity of the object.
    color : float
        Colour value (e.g., band1 - band2).
    color_bands : str
        String indicating which bands were used (e.g., 'g-r').
    band : str
        Band in which luminosity is measured (e.g., 'g').

    Returns
    -------
    float
        Estimated mass in solar masses.
    '''
    if color_bands == 'g-r':
        a = {'g': -0.499, 'r': -0.306, 'i': -0.222}
        b = {'g': 1.519, 'r': 1.097, 'i': 0.864}
    if color_bands == 'g-i':
        a = {'g': -0.379, 'r': -0.220, 'i': -0.152}
        b = {'g': 0.914, 'r': 0.661, 'i': 0.518}
    if color_bands == 'r-i':
        a = {'g': -0.106, 'r': -0.022, 'i': 0.006}
        b = {'g': 1.982, 'r': 1.431, 'i': 1.114}

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
        (g-i) colour value.
    abs_mag_i : float
        i-band absolute magnitude.

    Returns
    -------
    float
        Estimated mass in solar masses.
    '''
    return 10**(1.15 + 0.70 * color_g_i - 0.4 * abs_mag_i)


def stellmass2halo(stellar_mass):
    '''
    Converts stellar mass to halo mass using the relation from Huang+20:
        log10(M_stellar) [Msun] = 0.589 * (log10(M_halo) - 13.5) + 11.844

    Parameters
    ----------
    stellar_mass : float
        Stellar mass in solar masses.

    Returns
    -------
    float
        Halo (virial) mass in solar masses.
    '''
    halomass = 10**(((np.log10(stellar_mass) - 11.844) / 0.589) + 13.5)
    return halomass


# =============================================================================
# Function for dual axes plotting
# =============================================================================

def kpc_arcsec_axis(ax, kpc_ticks, arcsec_ticks, arcsec_lim, size_label, size_ticks, scale, DMpc, arcsec_label=True, kpc_label=True, z=None, spatial_label='R'):
    from matplotlib.ticker import AutoMinorLocator

    kpc_ticks = np.array(kpc_ticks)
    arcsec_ticks = np.array(arcsec_ticks)
    # Convert kpc ticks to arcsec positions.
    kpc_pos = Kpc2Arcsec(kpc_ticks, DMpc)
    
    # Only keep ticks within the specified arcsec limits.
    valid_range = (kpc_pos > arcsec_lim[0]) & (kpc_pos < arcsec_lim[1])
    kpc_pos = kpc_pos[valid_range]
    kpc_ticks = kpc_ticks[valid_range]
    ax.set_xticks(kpc_pos)
    ax.set_xticklabels([])
    if kpc_label:
        ax.set_xlabel(spatial_label + " [kpc]", size=size_label)
        ax.set_xticklabels(kpc_ticks)

    ax.set_xlim(arcsec_lim)
    
    # Process arcsec ticks similarly.
    valid_range = (arcsec_ticks > arcsec_lim[0]) & (arcsec_ticks < arcsec_lim[1])
    arcsec_ticks = arcsec_ticks[valid_range]

    ax_t = ax.twiny()
    ax_t.tick_params(axis="both", direction="out")
    ax_t.tick_params(axis='x')
    ax_t.set_xlim(arcsec_lim)
    ax_t.set_xticks(arcsec_ticks)
    ax_t.set_xticklabels([])
    ax_t.xaxis.set_minor_locator(AutoMinorLocator())

    if arcsec_label:
        ax_t.set_xlabel(spatial_label + " [arcsec]", size=size_label)
        ax_t.set_xticklabels(arcsec_ticks)


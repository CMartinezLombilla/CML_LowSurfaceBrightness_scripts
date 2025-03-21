#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: C. Martinez-Lombilla
"""

import numpy as np 


pixscale = 0.168 
zp = 27
z = 0.21007943

A = 0.06636 # r band
k_corr =  0.294962
dimming = 10*np.log10(1+z) 

sma0 = 2
step = 0.25 
nbins = 27



# =============================================================================
# Class to convert from counts to surface brightness [mag/arcsec^2]
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


# =============================================================================
# Function: extracts the  profile using the distances map
# =============================================================================

def dist_ellipse_prof(masked_data, dist_map, sma0, step, nbins, sky, zp, pixscale, A=0, dimming=0, k_corr=0):
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
    sky : float
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
        error_bin = sky/np.sqrt(float(np.size(data_bin)))
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
       
    return Table(vec, names=names, meta={'Counts and SB Profile; r in pix': 'table info'})
 
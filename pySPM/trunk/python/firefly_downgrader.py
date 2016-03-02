"""
Downgrades the resolution

"""
import numpy as np
from os import remove
import os.path
import astropy.io.fits as pyfits
import astropy
import astropy.constants
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline
from firefly_instrument import *
#-----------------------------------------------------------------------------


def downgrade(wave,flux,deltal_in,sigma_galaxy,wave_instrument,r_instrument):

    """
    Adapted from the manga DAP downgrader from Kyle Westfall.

    Downgrades an input spectrum to a given galaxy velocity dispersion
    using the input SEDs resolution and the resolution of the observation.

    Returns flux of downgraded SED.
    """

    
    # Define some constants
    # cnst = constants()

    sig2fwhm        = 2.0 * np.sqrt(2.0 * np.log(2.0))
    c               = 299792458.0 / 1000.0

    fwhm    = deltal_in/wave*c
    sigma   = fwhm/sig2fwhm
    sres    = wave/deltal_in


    new_sig     = np.zeros(wave.shape, dtype=np.float64)
    # match wavelength between model and instrument to downgrade
    def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return idx,array[idx]

    for wi,w in enumerate(wave):
        index,value = find_nearest(wave_instrument,w)
        sig_instrument = c/r_instrument[index]/sig2fwhm
        new_sig[wi] = np.sqrt(sigma_galaxy**2.0+sig_instrument**2.0)

    new_fwhm    = sig2fwhm*new_sig
    new_sres    = c/new_fwhm

    if len(wave)<5:
        raise ValueError("Not enough wavelength points...!")

    a_wave = wave[2]-wave[1]
    b_wave = wave[3]-wave[2]

    if b_wave-a_wave < 0.000001*a_wave:
        log_wave = False
    else:
        log_wave = True


    new_flux, matched_sres, sigma_offset, new_mask = match_spectral_resolution(wave, flux, sres, wave, new_sres, min_sig_pix=0.0, log10=log_wave, new_log10=log_wave)
    

    return new_flux

"""
.. class:: GalaxySurveyMUSE

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class GalaxySurveyMUSE is dedicated to handling MUSE survey and the class GalaxySpectrumMUSE to handling its spectra.

"""
from os.path import join
import os
import numpy as n
import astropy.io.fits as fits
from scipy.interpolate import interp1d
from MiscellanousFunctionsLibrary import *
import astropy.units as u

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p

class GalaxySurveyMUSE:
    """
    Loads the environement proper to the MUSE survey :
     * Defines all the proper paths in the database,
         * Opens the catalog 
     
    :param redshift_catalog: name of the MUSE redshift catalog (path to the fits file)
    :param calibration: if the class is loaded with intention of flux calibrating the MUSE data (boolean)"""
    def __init__(self,redshift_catalog="Catalog.spectra_MACS1931.fits"):
        self.redshift_catalog = redshift_catalog
        self.database_dir = os.environ['DATA_DIR']
        self.muse_dir = join(self.database_dir,"MUSE")
        self.muse_catalog_dir = join(self.muse_dir,"catalogs")
        self.muse_spectra_dir = join(self.muse_dir,"spectra")
        hd = fits.open(join(self.muse_catalog_dir,self.redshift_catalog))
        self.catalog = hd[1].data
        hd.close()

    def computeLineLuminosity(self,line,distanceCorrection):
        """ computes the line luminosities for the line list given.
        :param catalog: fits catalog containing redshift, EBV and line fluxes
        :param line:
        """
        ebvCorrection=n.array([ 10**(0.4 *self.catalog['SFD_EBV'][i] * CalzettiLaw((1 + self.catalog['ZBEST'][i]) * line[1])) for i in range(len(self.catalog['ZBEST']))])
        flux=ebvCorrection*self.catalog[line[2]+'_flux']*u.erg/u.cm**2/u.s
        Luminosity=fits.Column(name=line[2]+"_luminosity",format="D", unit="erg/s", array=distanceCorrection*flux )
        LuminosityErr=fits.Column(name=line[2]+"_luminosityErr",format="D", unit="erg/s", array= self.catalog[line[2]+'_fluxErr']/ self.catalog[line[2]+'_flux']* distanceCorrection *flux)
        return Luminosity, LuminosityErr



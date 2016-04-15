"""
.. class:: GalaxySpectrumDEEP2

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class GalaxySpectrumDEEP2 is dedicated to handling DEEP2 spectra

"""
from os.path import join
import os
import numpy as n
import astropy.io.fits as fits
import glob
from scipy.optimize import curve_fit
from GalaxySurveyMUSE import *

class GalaxySpectrumMUSE:
    """
    Loads the environement proper to the MUSE survey.

    Two modes of operation : flux calibration or line fitting
            
    :param catalog_entry: an entry of the MUSE catalog
    :param survey: survey python class
    :param calibration: if the class is loaded with intention of flux calibrating the MUSE data.
    :param lineFits: if the class is loaded with intention of fitting line fluxes on the MUSE spectra."""
    def __init__(self,catalog_entry, survey=GalaxySurveyMUSE()):
        self.catalog_entry=catalog_entry
        
        self.database_dir = os.environ['DATA_DIR']
        self.muse_dir = join(self.database_dir,"MUSE")
        self.muse_catalog_dir = join(self.muse_dir,"catalogs")
        self.muse_spectra_dir = join(self.muse_dir,"spectra")
        self.survey = survey
		self.path_to_spectrum = glob.glob(join(self.muse_spectra_dir , "spec_"+self.catalog_entry['SpecName']+".txt"))
        print "path to spectrum", self.path_to_spectrum

        
    def openObservedSpectrum(self):
        """Loads the observed spectrum in counts."""
        self.wavelength, self.fluxl, self.fluxlErr = n.loadtxt(self.path_to_spectrum[0])




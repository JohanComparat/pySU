"""
.. class:: GalaxySpectrumVVDS

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class GalaxySpectrumVVDS is dedicated to handling VVDS spectra

"""
from os.path import join
import os
import numpy as n
import astropy.io.fits as fits
import glob

class GalaxySpectrumVVDS:
        """
        Loads the environement proper to the vvds survey.

        Two modes of operation : flux calibration or line fitting
                
        :param catalog_entry: an entry of the vvds catalog
        :param calibration: if the class is loaded with intention of flux calibrating the vvds data.
        :param lineFits: if the class is loaded with intention of fitting line fluxes on the vvds spectra.
	"""
	def __init__(self,catalog_entry,lineFits=False):
		self.catalog_entry=catalog_entry
		self.database_dir = os.environ['DATA_DIR']
		self.vvds_dir = join(self.database_dir,"VVDS")
		self.vvds_catalog_dir = join(self.vvds_dir,"catalogs")
		self.vvds_spectra_dir = join(self.vvds_dir,"spectra")

	def openObservedSpectrum(self):
		"""
		reads a VVDS pectrum
		returns the wavelength, the flux and the error on the flux and two arrays for masking purpose
		"""
		spL=glob.glob(join(self.vvds_spectra_dir,"sc_*" + str(self.catalog_entry['NUM']) + "*atm_clean.fits"))
		#print spL
		if len(spL)==1 :
			specFileName=spL[0]
			spectraHDU=fits.open(specFileName)
			wl=spectraHDU[0].header['CRVAL1'] + spectraHDU[0].header['CDELT1'] * n.arange(2,spectraHDU[0].header['NAXIS1']+2)
			fl=spectraHDU[0].data[0]
			noiseFileName=glob.glob(join(self.vvds_spectra_dir+$,"sc_*"+str(self.catalog_entry['NUM'])+"*noise.fits"))[0]
			noiseHDU=fits.open(noiseFileName)
			flErr=noiseHDU[0].data[0]
			self.wavelength,self.fluxl,self.fluxlErr=wl,fl,flErr
		else :
			self.wavelength,self.fluxl,self.fluxlErr= [-1,-1.],[-1,-1.],[-1,-1.]




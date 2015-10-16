"""
.. class:: GalaxySpectrumVIPERS

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class GalaxySpectrumVIPERS is dedicated to handling VIPERS spectra

"""
from os.path import join
import os
import numpy as n
import astropy.io.fits as fits
import glob

class GalaxySpectrumVIPERS:
        """
        Loads the environement proper to the vipers survey.

        Two modes of operation : flux calibration or line fitting
                
        :param catalog_entry: an entry of the vipers catalog
        :param calibration: if the class is loaded with intention of flux calibrating the vipers data.
        :param lineFits: if the class is loaded with intention of fitting line fluxes on the vipers spectra.
	"""
	def __init__(self,catalog_entry,calibration=True,lineFits=False):
		self.catalog_entry=catalog_entry
		self.database_dir = os.environ['DATA_DIR']
		self.vipers_dir = join(self.database_dir,"VIPERS")
		self.vipers_catalog_dir = join(self.vipers_dir,"catalogs")
		self.vipers_spectra_dir = join(self.vipers_dir,"spectra")

	def openObservedSpectrum(self):
		"""
		reads a VIPERS pectrum
		returns the wavelength, the flux and the error on the flux and two arrays for masking purpose
		"""
		self.field='W'+self.catalog_entry['id_IAU'][7]
		specFileName=join(self.vipers_spectra_dir,"VIPERS_"+ self.field+ "_PDR1_SPECTRA_1D",self.catalog_entry['id_IAU'][:6]+"_"+self.catalog_entry['id_IAU'][7:]+".fits")
		spectraHDU=fits.open(specFileName)
		wlA=spectraHDU[1].data['WAVES']
		flA=spectraHDU[1].data['FLUXES']
		flErrA=spectraHDU[1].data['NOISE']
		edit=spectraHDU[1].data['EDIT']
		mask=spectraHDU[1].data['MASK']
		self.wavelength,self.fluxl,self.fluxlErr= wlA[(mask==0)], flA[(mask==0)], flErrA[(mask==0)]
		spectraHDU.close()





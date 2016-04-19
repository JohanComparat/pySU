"""
.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

*General purpose*:

The class GalaxySpectrumFIREFLY is dedicated to handling spectra to be fed to FIREFLY for fitting its stellar population

*Imports*::

	import numpy as np
	import astropy.io.fits as pyfits
	import glob
	from firefly_instrument import *
	from firefly_dust import *
	from firefly_fitter import *
	from firefly_library import *


"""
import numpy as np
import astropy.io.fits as pyfits
import glob

from firefly_instrument import *
from firefly_dust import *
from firefly_fitter import *
from firefly_library import *


class GalaxySpectrumFIREFLY:
	"""
	Loads the environnement to transform observed spectra into the input for FIREFLY. 
	
	Currently SDSS spectra, speclite format is handled as well as stacks from the VVDS and the DEEP2 galaxy surveys.

	:param path_to_spectrum: path to the spectrum
	:param hpf_mode: models the dust attenuation observed in the spectrum using high pass filter.
	:param milky_way_reddening: True if you want to correct from the Milky way redenning using the Schlegel 98 dust maps.
	:param stack_type: Optional. If you are fitting stack, choose the origins of the stacks : "DEEP2" or "VVDS"
	"""
	def __init__(self,path_to_spectrum, milky_way_reddening=True , hpf_mode = 'on', stack_resolution = 8686.):
		self.path_to_spectrum=path_to_spectrum
		self.milky_way_reddening = milky_way_reddening
		self.stack_resolution = stack_resolution 

	def openObservedMuseSpectrum(self, catalog):
		"""Loads the observed spectrum in counts."""
		self.wavelength, self.flux, self.error = n.loadtxt(self.path_to_spectrum, unpack=True)
		self.bad_flags = np.ones(len(self.wavelength))
		bad_data = np.isnan(self.flux) | np.isinf(self.flux) | (self.flux <= 0.0) | np.isnan(self.error) | np.isinf(self.error)
		# removes the bad data from the spectrum 
		self.flux[bad_data] 	= 0.0
		self.error[bad_data] 	= np.max(self.flux) * 99999999999.9
		self.bad_flags[bad_data] = 0

		self.r_instrument = np.zeros(len(self.wavelength))
		for wi,w in enumerate(self.wavelength):
			if w<6000:
				self.r_instrument[wi] = (2270.0-1560.0)/(6000.0-3700.0)*w + 420.0 
			else:
				self.r_instrument[wi] = (2650.0-1850.0)/(9000.0-6000.0)*w + 250.0 

		self.redshift = catalog['FINAL_Z']
		self.vdisp = 100 # catalog['VDISP']
		self.restframe_wavelength = self.wavelength / (1.0+self.redshift)

		self.trust_flag = 1
		self.objid = 0

		if self.milky_way_reddening :
			# gets the amount of MW reddening on the models
			self.ebv_mw = get_dust_radec(self.ra,self.dec,'ebv')
		else:
			self.ebv_mw = 0.0

	def openObservedSDSSSpectrum(self):
		"""
		It reads an SDSS spectrum and provides the input for the firefly fitting routine.
		:param path_to_spectrum:
		:param sdss_dir: directory with the observed spectra
		:param milky_way_reddening: True or False if you want to correct the redenning of the Milky way.
		:param hpf_mode: 'on' high pass filters the data to correct from dust in the galaxy.

		In this aims, it stores the following data in the object :
		* hdu list from the spec lite
		* SED data : wavelength (in angstrom), flux, error on the flux (in 10^{-17} erg/cm2/s/Angstrom, like the SDSS spectra)
		* Metadata :
			* ra : in degrees J2000
			* dec : in degrees J2000
			* redshift : best fit
			* vdisp : velocity dispersion in km/s
			* r_instrument : resolution of the instrument at each wavelength observed
			* trust_flag : 1 or True if trusted 
			* bad_flags : ones as long as the wavelength array, filters the pixels with bad data
			* objid : object id optional : set to 0
		"""
		self.hdulist = pyfits.open(self.path_to_spectrum)
		self.ra = self.hdulist[0].header['RA']
		self.dec = self.hdulist[0].header['DEC']

		self.wavelength = 10**self.hdulist[1].data['loglam']
		self.flux = self.hdulist[1].data['flux']
		self.error = self.hdulist[1].data['ivar']**(-0.5)
		self.bad_flags = np.ones(len(self.wavelength))
		bad_data = np.isnan(self.flux) | np.isinf(self.flux) | (self.flux <= 0.0) | np.isnan(self.error) | np.isinf(self.error)
		# removes the bad data from the spectrum 
		self.flux[bad_data] 	= 0.0
		self.error[bad_data] 	= np.max(self.flux) * 99999999999.9
		self.bad_flags[bad_data] = 0

		self.r_instrument = np.zeros(len(self.wavelength))
		for wi,w in enumerate(self.wavelength):
			if w<6000:
				self.r_instrument[wi] = (2270.0-1560.0)/(6000.0-3700.0)*w + 420.0 
			else:
				self.r_instrument[wi] = (2650.0-1850.0)/(9000.0-6000.0)*w + 250.0 

		self.redshift = self.hdulist[2].data['Z'][0]
		self.vdisp = self.hdulist[2].data['VDISP'][0]
		self.restframe_wavelength = self.wavelength / (1.0+self.redshift)

		self.trust_flag = 1
		self.objid = 0

		if self.milky_way_reddening :
			# gets the amount of MW reddening on the models
			self.ebv_mw = get_dust_radec(self.ra,self.dec,'ebv')
		else:
			self.ebv_mw = 0.0


	def openObservedStack(self):
		"""
		It reads an Stack spectrum from the LF analysis and provides the input for the firefly fitting routine.
		:param path_to_spectrum:
		:param sdss_dir: directory with the observed spectra
		:param milky_way_reddening: True or False if you want to correct the redenning of the Milky way.
		:param hpf_mode: 'on' high pass filters the data to correct from dust in the galaxy.

		In this aims, it stores the following data in the object :
		* hdu list from the spec lite
		* SED data : wavelength (in angstrom), flux, error on the flux (in 10^{-17} erg/cm2/s/Angstrom, like the SDSS spectra)
		* Metadata :
			* ra : in degrees J2000
			* dec : in degrees J2000
			* redshift : best fit
			* vdisp : velocity dispersion in km/s
			* r_instrument : resolution of the instrument at each wavelength observed
			* trust_flag : 1 or True if trusted 
			* bad_flags : ones as long as the wavelength array, filters the pixels with bad data
			* objid : object id optional : set to 0
		"""
		self.hdulist = pyfits.open(self.path_to_spectrum)
		self.ra = 0. #self.hdulist[0].header['RA']
		self.dec = 0. #self.hdulist[0].header['DEC']
		self.redshift = float(self.path_to_spectrum.split('-')[-1].split('_')[0][1:])
		self.restframe_wavelength = self.hdulist[1].data['wavelength']

		self.wavelength = self.restframe_wavelength * (1. + self.redshift)
		self.flux = self.hdulist[1].data['meanWeightedStack'] * 1e17
		self.error = self.hdulist[1].data['jackknifStackErrors'] * 1e17
		self.bad_flags = np.ones(len(self.restframe_wavelength))
		Nstacked = float(self.path_to_spectrum.split('-')[-1].split('_')[3])

		N_angstrom_masked = 20
		lines_mask = ((self.restframe_wavelength > 3738 - N_angstrom_masked) & (self.restframe_wavelength < 3738 + N_angstrom_masked)) | ((self.restframe_wavelength > 5007 - N_angstrom_masked) & (self.restframe_wavelength < 5007 + N_angstrom_masked)) | ((self.restframe_wavelength > 4861 - N_angstrom_masked) & (self.restframe_wavelength < 4861 + N_angstrom_masked)) | ((self.restframe_wavelength > 6564 - N_angstrom_masked) & (self.restframe_wavelength < 6564 + N_angstrom_masked)) | ( self.hdulist[1].data['NspectraPerPixel'] < Nstacked * 0.8 ) | (self.flux==-9999.99)

		self.restframe_wavelength = self.restframe_wavelength[(lines_mask==False)] 
		self.wavelength = self.wavelength[(lines_mask==False)] 
		self.flux = self.flux[(lines_mask==False)] 
		self.error = self.error[(lines_mask==False)] 
		self.bad_flags = self.bad_flags[(lines_mask==False)] 

		bad_data 	= np.isnan(self.flux) | np.isinf(self.flux) | (self.flux <= 0.0) | np.isnan(self.error) | np.isinf(self.error) 
		# removes the bad data from the spectrum 
		self.flux[bad_data] 	= 0.0
		self.error[bad_data] 	= np.max(self.flux) * 99999999999.9
		self.bad_flags[bad_data] = 0

		self.r_instrument = np.ones_like(np.arange(len(self.restframe_wavelength))) * self.stack_resolution * 2.

		self.vdisp = 70. # km/s

		self.trust_flag = 1
		self.objid = 0

		if self.milky_way_reddening :
			# gets the amount of MW reddening on the models
			self.ebv_mw = get_dust_radec(self.ra,self.dec,'ebv')
		else:
			self.ebv_mw = 0.0

		print"there are", len(self.wavelength),"data points at redshift",self.redshift," between", np.min(self.wavelength[bad_data==False]), np.max(self.wavelength[bad_data==False]), "Angstrom.", np.min(self.restframe_wavelength[bad_data==False]), np.max(self.restframe_wavelength[bad_data==False]), "Angstrom in the rest frame."


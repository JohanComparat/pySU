"""
.. class:: GalaxySurveyVVDS

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class GalaxySurveyVVDS is dedicated to handling VVDS survey.

"""
import numpy as n
import astropy.io.fits as fits
from MiscellanousFunctionsLibrary import *
import astropy.units as u

class GalaxySurveyVVDS:
	"""
	Loads the environement proper to the VVDS survey :
	 * Defines all the proper paths in the database,
         * Opens the catalog,

        :param redshift_catalog: name of the VVDS redshift catalog (path to the fits file)
	"""
	def __init__(self,redshift_catalog="VVDS_WIDE_summary.fits"):
		self.redshift_catalog = redshift_catalog
		self.database_dir = "/home/comparat/database/"
		self.vvds_dir = self.database_dir+"VVDS/"
		self.vvds_catalog_dir = self.vvds_dir+"catalogs/"
		self.vvds_spectra_dir = self.vvds_dir+"spectra/"
		hd = fits.open(self.vvds_catalog_dir+self.redshift_catalog)
		self.catalog = hd[1].data
		hd.close()

	def computeLineLuminosity(self,line,distanceCorrection):
		""" computes the line luminosities for the line list given.
		:param catalog: fits catalog containing redshift, EBV and line fluxes
		:param line:
		"""
		ebvCorrection=n.array([ 10**(0.4 *self.catalog['EBV_MW'][i] * CalzettiLaw((1 + self.catalog['Z'][i]) * line[1])) for i in range(len(self.catalog['Z']))])
		correctionAperture = 1. / self.catalog['fo']
		flux=ebvCorrection* correctionAperture * self.catalog[line[2] +'_flux']* u.erg/ u.cm**2 /u.s 
		Luminosity=fits.Column(name=line[2]+"_luminosity",format="D", unit="erg/s", array=distanceCorrection*flux )
		LuminosityErr=fits.Column(name=line[2]+"_luminosityErr",format="D", unit="erg/s", array= self.catalog[line[2]+'_fluxErr']/ self.catalog[line[2]+'_flux']* distanceCorrection *flux)
		return Luminosity, LuminosityErr



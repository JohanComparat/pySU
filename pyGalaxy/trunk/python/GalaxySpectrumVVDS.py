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
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p
from LineFittingLibrary import *
lfl = LineFittingLibrary()
from filterList import *

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
			noiseFileName=glob.glob(join(self.vvds_spectra_dir,"sc_*"+str(self.catalog_entry['NUM'])+"*noise.fits"))[0]
			noiseHDU=fits.open(noiseFileName)
			flErr=noiseHDU[0].data[0]
			self.wavelength,self.fluxl,self.fluxlErr=wl,fl,flErr
		else :
			self.wavelength,self.fluxl,self.fluxlErr= [-1,-1.],[-1,-1.],[-1,-1.]

	def plotFit(self, outputFigureNameRoot):
		"""
		Plots the spectrum and the line fits in a few figures
		"""
		ifl = lfl.flambda(self.catalog_entry['MAGI'], lambIcfht)
		ifl_max = lfl.flambda(self.catalog_entry['MAGI']+self.catalog_entry['MAGERR_AUTO_I_1'], lambIcfht)
		ifl_min = lfl.flambda(self.catalog_entry['MAGI']-self.catalog_entry['MAGERR_AUTO_I_1'], lambIcfht)

		rfl = lfl.flambda(self.catalog_entry['MAG_R_CFHTLS'], lambRcfht)
		rfl_max = lfl.flambda(self.catalog_entry['MAG_R_CFHTLS']+self.catalog_entry['MAGERR_AUTO_R_1'], lambRcfht)
		rfl_min = lfl.flambda(self.catalog_entry['MAG_R_CFHTLS']-self.catalog_entry['MAGERR_AUTO_R_1'], lambRcfht)
		
		p.figure(1,(12,4))
		p.axes([0.1,0.2,0.85,0.75])
		p.errorbar(self.wavelength,self.fluxl,yerr = self.fluxlErr, linewidth=1, label='spectrum')
		p.plot([lambIcfht,lambIcfht,lambIcfht],[ifl_min,ifl,ifl_max], 'b,', label = 'magnitudes')
		p.plot([lambRcfht, lambRcfht, lambRcfht], [rfl_min, rfl, rfl_max], 'b,')
		p.xlabel('wavelength [A]')
		p.ylabel(r'f$_\lambda$ [erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
		p.yscale('log')
		p.ylim((1e-20, 1e-18))
		p.legend(fontsize=12) 
		p.savefig( outputFigureNameRoot + "-all.png" )
		p.clf()

		a0 = self.catalog_entry['O2_3728_a0']
		continu= self.catalog_entry['O2_3728_continu']
		aas =n.arange(self.catalog_entry['O2_3728_a0']-30, self.catalog_entry['O2_3728_a0']+30,0.1)
		flMod=lambda aa,sigma,F0,sh :continu+ lfl.gaussianLineNC(aa,sigma,(1-sh)*F0,a0)+lfl.gaussianLineNC(aa,sigma,sh*F0,a0)
		model = flMod(aas, self.catalog_entry['O2_3728_sigma'], self.catalog_entry['O2_3728_flux'],0.58 )# self.catalog_entry['O2_3728_share'])
		
		p.figure(2,(4,4))
		p.axes([0.1,0.2,0.85,0.75])
		p.errorbar(self.wavelength,self.fluxl,yerr = self.fluxlErr)
		p.axvline(self.catalog_entry['O2_3728_a0'],color='k', ls='dashed', label= 'obs')
		p.plot(aas, model,'g',label='model')
		p.xlabel('wavelength [A]')
		p.ylabel(r'f$_\lambda$ [erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
		p.yscale('log')
		p.ylim((1e-20, 1e-18))
		p.xlim(( self.catalog_entry['O2_3728_a0']-50, self.catalog_entry['O2_3728_a0']+50))
		p.legend(fontsize=12, loc=4) 
		p.savefig( outputFigureNameRoot + "-O2_3728.png")
		p.clf()

		a0 = self.catalog_entry['O3_5007_a0']
		continu= self.catalog_entry['O3_5007_continu']
		aas =n.arange(self.catalog_entry['O3_5007_a0']-30, self.catalog_entry['O3_5007_a0']+30,0.1)
		flMod=lambda aa,sigma,F0: lfl.gaussianLine(aa,sigma,F0,a0,continu)
		model = flMod(aas, self.catalog_entry['O3_5007_sigma'], self.catalog_entry['O3_5007_flux'])
		
		p.figure(2,(4,4))
		p.axes([0.1,0.2,0.85,0.75])
		p.errorbar(self.wavelength,self.fluxl,yerr = self.fluxlErr)
		p.axvline(self.catalog_entry['O3_5007_a0'],color='k', ls='dashed', label= 'obs')
		p.plot(aas, model,'g',label='model')
		p.xlabel('wavelength [A]')
		p.ylabel(r'f$_\lambda$ [erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
		p.yscale('log')
		p.ylim((1e-20, 1e-18))
		p.xlim(( self.catalog_entry['O3_5007_a0']-50, self.catalog_entry['O3_5007_a0']+50))
		p.legend(fontsize=12, loc=4) 
		p.savefig( outputFigureNameRoot + "-O3_5007.png")
		p.clf()



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
from lineListAir import *

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

	def plotFit(self, outputFigureNameRoot, ymin = 1e-19, ymax = 1e-17):
		"""
		Plots the spectrum and the line fits in a few figures
		"""
		ifl = lfl.flambda(self.catalog_entry['MAGI'], lambIcfht)
		ifl_max = lfl.flambda(self.catalog_entry['MAGI']+self.catalog_entry['MAGERR_AUTO_I_1'], lambIcfht)
		ifl_min = lfl.flambda(self.catalog_entry['MAGI']-self.catalog_entry['MAGERR_AUTO_I_1'], lambIcfht)

		rfl = lfl.flambda(self.catalog_entry['MAG_R_CFHTLS'], lambRcfht)
		rfl_max = lfl.flambda(self.catalog_entry['MAG_R_CFHTLS']+self.catalog_entry['MAGERR_AUTO_R_1'], lambRcfht)
		rfl_min = lfl.flambda(self.catalog_entry['MAG_R_CFHTLS']-self.catalog_entry['MAGERR_AUTO_R_1'], lambRcfht)
		
		ok = (self.fluxl >0 ) & (self.fluxl > 2* self.fluxlErr)
		p.figure(1,(12,4))
		p.axes([0.1,0.2,0.85,0.75])
		p.errorbar(self.wavelength[ok],self.fluxl[ok],yerr = self.fluxlErr[ok], linewidth=1, alpha= 0.4, label='spectrum')
		p.plot([lambIcfht,lambIcfht,lambIcfht],[ifl_min,ifl,ifl_max], 'r', label = 'magnitudes', lw=2)
		p.plot([lambRcfht, lambRcfht, lambRcfht], [rfl_min, rfl, rfl_max], 'r', lw=2)
		p.xlabel('wavelength [A]')
		p.ylabel(r'f$_\lambda$ [erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
		p.yscale('log')
		p.ylim((ymin, ymax))
		gl = p.legend(loc=0,fontsize=12)
		gl.set_frame_on(False)
		p.savefig( outputFigureNameRoot + "-all.png" )
		p.clf()

		a0_1 = (1+self.catalog_entry['Z'])*O2_3727
		a0_2 = (1+self.catalog_entry['Z'])*O2_3729
		continu= self.catalog_entry['O2_3728_continu']
		aas =n.arange(self.catalog_entry['O2_3728_a0']-70, self.catalog_entry['O2_3728_a0']+70,0.1)
		flMod=lambda aa,sigma,F0,sh :continu+ lfl.gaussianLineNC(aa,sigma,(1-sh)*F0,a0_1)+lfl.gaussianLineNC(aa,sigma,sh*F0,a0_2)
		model = flMod(aas, self.catalog_entry['O2_3728_sigma'], self.catalog_entry['O2_3728_flux'],0.58 )# self.catalog_entry['O2_3728_share'])
		
		p.figure(2,(4,4))
		p.axes([0.21,0.2,0.78,0.7])
		p.errorbar(self.wavelength,self.fluxl,yerr = self.fluxlErr)
		p.plot(aas, model,'g',label='model', lw=2)
		p.xlabel('wavelength [A]')
		p.ylabel(r'f$_\lambda$ [erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
		p.yscale('log')
		p.ylim((ymin, ymax))
		p.xlim(( self.catalog_entry['O2_3728_a0']-100, self.catalog_entry['O2_3728_a0']+100))
		gl = p.legend(loc=0,fontsize=12)
		gl.set_frame_on(False)
		p.title('[OII] 3727')
		p.savefig( outputFigureNameRoot + "-O2_3728.png")
		p.clf()

		a0 = self.catalog_entry['O3_5007_a0']
		continu= self.catalog_entry['O3_5007_continu']
		aas =n.arange(self.catalog_entry['O3_5007_a0']-70, self.catalog_entry['O3_5007_a0']+70,0.1)
		flMod=lambda aa,sigma,F0: lfl.gaussianLine(aa,sigma,F0,a0,continu)
		model = flMod(aas, self.catalog_entry['O3_5007_sigma'], self.catalog_entry['O3_5007_flux'])
		
		p.figure(2,(4,4))
		p.axes([0.21,0.2,0.78,0.7])
		p.errorbar(self.wavelength,self.fluxl,yerr = self.fluxlErr)
		p.plot(aas, model,'g',label='model', lw =2)
		p.xlabel('wavelength [A]')
		p.ylabel(r'f$_\lambda$ [erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
		p.yscale('log')
		p.ylim((ymin, ymax))
		p.xlim(( self.catalog_entry['O3_5007_a0']-100, self.catalog_entry['O3_5007_a0']+100))
		gl = p.legend(loc=0,fontsize=12)
		gl.set_frame_on(False)
		p.title('[OIII] 5007')
		p.savefig( outputFigureNameRoot + "-O3_5007.png")
		p.clf()



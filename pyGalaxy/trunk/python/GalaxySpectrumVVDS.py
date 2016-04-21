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
from lineFittingLibrary import *
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
		ifl = flambda(self.catalog_entry['MAGI'], lambIcfht)
		rfl = flambda(self.catalog_entry['MAG_R_CFHTLS'], lambRcfht)
		p.figure(1,(12,4))
		p.axes([0.1,0.2,0.85,0.75])
		p.errorbar(self.wavelength,self.fluxl,yerr = self.fluxlErr, linewidth=1)
		p.plot()
		p.xlabel('wavelength [A]')
		p.ylabel(r'f$_\lambda$ [erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
		p.yscale('log')
		p.legend(fontsize=12) 
		p.savefig( outputFigureNameRoot + "-all.png" )
		p.clf()

		p.figure(2,(12,4))
		p.axes([0.1,0.2,0.85,0.75])
		p.errorbar(self.wavelength,self.fluxl,yerr = self.fluxlErr,label= "z="+str(redshift)+", "+str(Lmin) +"< log(L"+line+")<"+ str(Lmax),linewidth=1)
		p.axvline(4861,color='k', ls='dashed')
		p.axvline(4341,color='k', ls='dashed')
		try :
			aa=str(n.round(hdus[0].header['H1_4862_flux_nc']/hdus[0].header['H1_4341_flux_nc'] , 3))
			bb=str(n.round(hdus[0].header['EBV_4862_4341'],3))
			p.text(4331+50,1e-17,r"GP. H$\beta/\gamma$="+aa+", EBV=" +bb)
		except KeyError:
			pass

		p.xlabel('wavelength [A]')
		p.ylabel(r'f$_\lambda$ [erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
		p.yscale('log')
		p.xlim(( 4331, 4871 ))
		p.legend(fontsize=12, loc=4) 
		p.savefig( join( os.environ['SPECTRASTACKS_DIR'], "plots", "models", modeledStackFile.split('/')[-1] + "-H1lineDecr.png"))
		p.clf()



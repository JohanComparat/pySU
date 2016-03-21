"""
.. class:: SpectraStacking

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class SpectraStacking is dedicated to stacking spectra

"""
import os 
#import astropy.cosmology as co
#cosmo=co.Planck15 # co.FlatLambdaCDM(H0=70,Om0=0.3)
import astropy.io.fits as fits
import numpy as n
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.stats import scoreatpercentile

from MiscellanousFunctionsLibrary import *
from HandleReducedELGPlate import *
#SpectraStacking("/home/comparat/database/DEEP2/products/emissionLineLuminosityFunctions/O2_3728/O2_3728-DEEP2-z0.925.fits")

class SpectraStacking:
	"""
	The model luminosity function class
	:param LF_file: fits file generated with a LF.
	:param Resolution: Resolution
	:param Nspec: Initial guess of the parameters
	:param outputDirectory: where to output the fits
	:param fileName: file name where things are saved.
	"""
	def __init__(self, LF_file,Nspec= 400, dLambda = 0.0001, dV=-9999.99):
		self.LF_file = LF_file
		self.dLambda = dLambda
		self.wave= 10**n.arange(3.1760912590556813, 4.0211892990699383, dLambda) # 1500,10500
		self.R = int(1/n.mean((self.wave[1:] -self.wave[:-1])/ self.wave[1:]))
		self.dV = dV
		self.Nspec = Nspec
		self.catalog_entries=fits.open(self.LF_file)[1].data
		print len(self.catalog_entries)

	def stack_function(self,specMatrix,specMatrixWeight):
		"""Creates the stack.
		:param specMatrix: matrix of observed spectra
		:param specMatrixWeight: matrix of the statistical weights used in the LF.
		"""
		stackMed = n.ones_like(n.empty(len(self.wave)))*self.dV
		stackMean = n.ones_like(n.empty(len(self.wave)))*self.dV
		stackMeanWeighted = n.ones_like(n.empty(len(self.wave)))*self.dV
		stackVar = n.ones_like(n.empty(len(self.wave)))*self.dV
		stackN = n.ones_like(n.empty(len(self.wave)))*self.dV
		jackknifes = n.ones_like(n.empty((len(self.wave),10)))*self.dV
		for i in range(len(specMatrix.T)):
				pt=specMatrix.T[i]
				wt=specMatrixWeight.T[i]
				sel=(pt!=self.dV)
				# jackknife sub-sampling
				rd=n.random.random(len(pt))
				aim=n.arange(0,1.01,0.1)
				jks=n.array([ (rd>aim[jj])&(rd<aim[jj+1]) for jj in range(len(aim)-1) ])
				if len(pt[sel])>1:
						stackMed[i] = n.median(pt[sel])
						stackMean[i] = n.mean(pt[sel])
						stackMeanWeighted[i] = n.average(pt[sel],weights=wt[sel])
						stackN[i] = len(pt[sel])
						inter = n.array([ n.median( pt[sel & (seK==False)] ) for seK in jks ])
						jackknifes[i] = inter
						stackVar[i] = n.std(inter)

		wavelength = fits.Column(name="wavelength",format="D", unit="Angstrom", array= self.wave)
		medianStack=fits.Column(name="medianStack",format="D", unit="erg/s/cm2/Angstrom", array= n.array(stackMed))
		meanStack=fits.Column(name="meanStack",format="D", unit="erg/s/cm2/Angstrom", array= n.array(stackMean))
		meanWeightedStack=fits.Column(name="meanWeightedStack",format="D", unit= "erg/s/cm2/Angstrom", array= n.array(stackMeanWeighted))
		jackknifeSpectra=fits.Column(name="jackknifeSpectra",format="10D", unit="erg/s/cm2/Angstrom", array= n.array(jackknifes))
		jackknifStackErrors=fits.Column(name="jackknifStackErrors",format="D", unit="erg/s/cm2/Angstrom", array= n.array(stackVar))
		NspectraPerPixel=fits.Column(name="NspectraPerPixel",format="D", unit="", array= n.array(stackN))
		return  wavelength, medianStack, meanStack, meanWeightedStack, jackknifStackErrors, jackknifeSpectra, NspectraPerPixel

	def convertSpectrum(self,redshift):
		"""
		Shifts the spectrum in the rest-frame and creates a spectrum with the sampling desired.
		:param redshift: redshift of the spectrum
		"""	
		nwave=self.wavelength/(1+redshift)

		inL=(self.wave>nwave.min())&(self.wave<nwave.max())
		outL=(inL==False)

		points=interp1d(nwave,nwave * self.fluxl)
		pts=points(self.wave[inL]) / self.wave[inL]
		res=n.ones_like(self.wave)*self.dV
		res[inL]=pts

		pointsErr=interp1d(nwave,nwave * self.fluxlErr)
		ptsErr=pointsErr(self.wave[inL]) / self.wave[inL]
		resErr=n.ones_like(self.wave)*self.dV
		resErr[inL]=ptsErr

		return res, resErr

	def stackEbossPlateSpectra(self,outPutFileName, g_min = 22,g_max=22.8, gr_min=-0.5, gr_max=1.5, rz_min=-0.5, rz_max = 1.5):
		"""
		Function that constructs the stacks for a listof spectra
		"""
		# loop over the file with N sorted with luminosity
		indexes = n.arange(len(self.catalog_entries['PLATE']))
		jumps = n.arange(0, len(self.catalog_entries['PLATE']), self.Nspec)
		for ii in range(len(jumps)-1):
			ids = indexes[jumps[ii]:jumps[ii+1]]
			specMatrix,specMatrixErr,specMatrixWeight=[],[],[]
			count=0
			Zdistrib = scoreatpercentile( self.catalog_entries[ids]['Z' ] , [0,25,50,75,100])
			print "stacks ",len(self.catalog_entries[ids]), "galaxies from eBOSS plates with redshifts (min, Q25, median, Q75, max)", Zdistrib
			for entry in self.catalog_entries[ids] :
				ObsPlate = HandleReducedELGPlate(entry['PLATE'],entry['MJD'])
				ObsPlate.loadPlate()
				print "Z = ", ObsPlate.zstruc['Z'][entry['FIBERID']-1], entry['Z']
				self.wavelength,self.fluxl,self.fluxlErr = ObsPlate.wavelength, ObsPlate.flux[entry['FIBERID']-1], ObsPlate.fluxErr[entry['FIBERID']-1]
				pts,ptsErr = self.convertSpectrum(entry['Z'])
				specMatrix.append(pts)
				specMatrixErr.append(ptsErr)
				weight=1.
				specMatrixWeight.append(n.ones_like(pts)*weight)
				count+=1

			specMatrixWeight=n.array(specMatrixWeight)
			specMatrix=n.array(specMatrix)
			specMatrixErr=n.array(specMatrixErr)
			print "now stacks"
			wavelength, medianStack, meanStack, meanWeightedStack, jackknifStackErrors, jackknifeSpectra, NspectraPerPixel = self.stack_function( specMatrix ,specMatrixWeight)
			cols = fits.ColDefs([wavelength, medianStack, meanStack, meanWeightedStack, jackknifStackErrors, jackknifeSpectra, NspectraPerPixel])
			tbhdu = fits.BinTableHDU.from_columns(cols)
			prihdr = fits.Header()
			prihdr['g_min'] = g_min
			prihdr['g_max'] = g_max
			prihdr['gr_min'] = gr_min
			prihdr['gr_max'] = gr_max
			prihdr['rz_min'] = rz_min
			prihdr['zr_max'] = rz_max
			prihdr['Z_min'] = Zdistrib[0]
			prihdr['Z_mean'] = Zdistrib[2]
			prihdr['Z_max'] = Zdistrib[-1]
			prihdu = fits.PrimaryHDU(header=prihdr)
			thdulist = fits.HDUList([prihdu, tbhdu])
			os.system('rm '+outPutFileName)
			thdulist.writeto(outPutFileName)


	def stackSdssMainSpectra(self, outPutFileName):
		"""
		Function that constructs the stacks for a listof spectra
		"""
		# loop over the file with N sorted with luminosity
		specDirSdssMain = "/uufs/chpc.utah.edu/common/home/sdss00/sdsswork/sdss/spectro/redux/26/spectra"
		indexes = n.argsort(self.catalog_entries['REDSHIFT'])
		jumps = n.arange(0, len(self.catalog_entries['PLATE']), self.Nspec)
		for ii in range(len(jumps)-1):
			ids = indexes[jumps[ii]:jumps[ii+1]]
			specMatrix,specMatrixErr,specMatrixWeight=[],[],[]
			count=0
			Zdistrib = scoreatpercentile( self.catalog_entries[ids]['REDSHIFT'] , [0,25,50,75,100])
			print "stacks ",len(self.catalog_entries[ids]), "galaxies from SDSS MAIN with redshifts (min, Q25, median, Q75, max)", Zdistrib
			for entry in self.catalog_entries[ids] :
				plate  = str(entry['PLATE']).zfill(4)
				specFile = join(specDirSdssMain, str(plate), "spec-"+str(plate) +"-" + str(entry['MJD'])+"-"+str(entry['FIBERID']).zfill(4)+".fits")
				print specFile
				hdulist = fits.open(specFile)
				sel = (hdulist[1].data['and_mask']==0)
				self.wavelength = 10.**hdulist[1].data['loglam'][sel]
				self.fluxl = hdulist[1].data['flux'][sel]
				self.fluxlErr = hdulist[1].data['ivar'][sel]**(-0.5)
				hdulist.close()
				pts,ptsErr = self.convertSpectrum(entry['REDSHIFT'])
				specMatrix.append(pts)
				specMatrixErr.append(ptsErr)
				weight=entry['WEIGHT_COMPLETENESS']*entry['WEIGHT_FC']
				specMatrixWeight.append(n.ones_like(pts)*weight)
				count+=1

			specMatrixWeight=n.array(specMatrixWeight)
			specMatrix=n.array(specMatrix)
			specMatrixErr=n.array(specMatrixErr)
			print "now stacks"
			wavelength, medianStack, meanStack, meanWeightedStack, jackknifStackErrors, jackknifeSpectra, NspectraPerPixel = self.stack_function( specMatrix ,specMatrixWeight)
			cols = fits.ColDefs([wavelength, medianStack, meanStack, meanWeightedStack, jackknifStackErrors, jackknifeSpectra, NspectraPerPixel])
			tbhdu = fits.BinTableHDU.from_columns(cols)
			prihdr = fits.Header()
			prihdr['Z_min'] = Zdistrib[0]
			prihdr['Z_mean'] = Zdistrib[2]
			prihdr['Z_max'] = Zdistrib[-1]
			prihdu = fits.PrimaryHDU(header=prihdr)
			thdulist = fits.HDUList([prihdu, tbhdu])
			os.system('rm '+outPutFileName)
			thdulist.writeto(outPutFileName)

"""
.. class:: SpectraStackingEBOSS

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class SpectraStacking is dedicated to stacking spectra from SDSS-IV eBOSS

"""
import os 
import astropy.io.fits as fits
import numpy as n
from scipy.interpolate import interp1d

class SpectraStackingEBOSS:
	"""
	The model luminosity function class
	:param in_file: file containing spectraidsto be stacked
	:param Resolution: Resolution
	:param out_file: where to output stacks
	"""
	def __init__(self, in_file, out_file, dLambda = 0.0001, dV=-9999.99, run2d = "v5_10_2"):
		print "input list:", in_file
		self.in_file = in_file
		self.plates, self.mjds, self.fiberids = n.loadtxt(self.in_file, unpack=True, dtype='str')
		self.out_file = out_file
		self.dLambda = dLambda
		self.wave= 10**n.arange(3.1760912590556813, 4.0211892990699383, dLambda) # 1500,10500
		self.R = int(1/n.mean((self.wave[1:] -self.wave[:-1])/ self.wave[1:]))
		print "R=", n.median(self.R)
		self.dV = dV
		self.survey="eBOSS"
		self.run2d = run2d
		self.run1d = self.run2d
		self.topdirBOSS = os.path.join(os.environ['BOSS_SPECTRO_REDUX'], run2d)

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

	def getSpectra(self, plate, mjd, fiber):
		# plate, mjd, fiber must be string !
		spfile = os.path.join(self.topdirBOSS, plate , 'spPlate-' + plate + '-' + mjd + '.fits')
		hdulist = fits.open(spfile)
		c0 = hdulist[0].header['coeff0']
		c1 = hdulist[0].header['coeff1']
		npix = hdulist[0].header['naxis1']
		self.wavelength = 10.**(c0 + c1 * n.arange(npix))  # here are the wavelengths
		flux = hdulist[0].data                                      # flux
		ivar = hdulist[1].data                                      # inverse variance
		hdulist.close()
		spZBfile = os.path.join(self.topdirBOSS, plate, self.run1d, 'spZbest-' + plate + '-' + mjd + '.fits')
		hdulist2 = fits.open(spZBfile)
		# now extract the spectrum
		i=int(fiber)-1
		redshift = hdulist2[1].data['Z'][i]
		hdulist2.close()
		self.fluxl = flux[i]
		iv=ivar[i]
		self.fluxlErr = iv**(-0.5)
		return redshift

	def stackSpectraSDSS(self):
		"""
		Function that constructs the stacks for a luminosity function. 
		It loops over the list of spectra given in the catalog of the LF. 
		First it sorts the catalog by the line luminosity. 
		And then stacks the first Nspec, then the next Nspec together.
		"""
		# loop over the file with N sorted with luminosity
		specMatrix, specMatrixErr, specMatrixWeight=[],[],[]
		
		for plate, mjd, fiber in zip(self.plates, self.mjds, self.fiberids):
				redshift = self.getSpectra(plate, mjd, fiber)
				pts,ptsErr = self.convertSpectrum(redshift)
				specMatrix.append(pts)
				specMatrixErr.append(ptsErr)
				weight=1.
				specMatrixWeight.append(n.ones_like(pts)*weight)

		specMatrixWeight=n.array(specMatrixWeight)
		specMatrix=n.array(specMatrix)
		specMatrixErr=n.array(specMatrixErr)
		print "now stacks"
		wavelength, medianStack, meanStack, meanWeightedStack, jackknifStackErrors, jackknifeSpectra, NspectraPerPixel = self.stack_function( specMatrix ,specMatrixWeight)
		cols = fits.ColDefs([wavelength, medianStack, meanStack, meanWeightedStack, jackknifStackErrors, jackknifeSpectra, NspectraPerPixel])
		tbhdu = fits.BinTableHDU.from_columns(cols)
		prihdr = fits.Header()
		prihdr['author'] = "JC"
		prihdr['survey'] = "eBOSS"
		prihdr['in_file'] = os.path.basename(self.in_file)[:-4]
		prihdr['Nspec'] = len(self.plates)
		prihdu = fits.PrimaryHDU(header=prihdr)
		thdulist = fits.HDUList([prihdu, tbhdu])
		if os.path.isfile(self.out_file):
			os.remove(self.out_file)
		print "stack written to", self.out_file
		thdulist.writeto(self.out_file)

"""
.. class:: SpectraStackingEBOSS

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class SpectraStacking is dedicated to stacking spectra from SDSS-IV eBOSS

Current version

"""
import os 
import astropy.io.fits as fits
import numpy as n
from scipy.interpolate import interp1d

maskLambda = n.loadtxt(os.path.join(os.environ['GIT_ARCHETYPES'],'data',"dr12-sky-mask.txt"), unpack=True)

get_path_to_spectrum_v5_10_0 = lambda plate, mjd, fiberid : os.path.join(os.environ['HOME'], 'SDSS', 'v5_10_0', 'spectra', str(int(plate)).zfill(4), "spec-"+str(int(plate)).zfill(4)+"-"+str(int(mjd)).zfill(5)+"-"+str(int(fiberid)).zfill(4)+".fits" )
get_path_to_spectrum_v5_10_7 = lambda plate, mjd, fiberid : os.path.join(os.environ['HOME'], 'SDSS', 'v5_10_7', 'spectra', str(int(plate)).zfill(4), "spec-"+str(int(plate)).zfill(4)+"-"+str(int(mjd)).zfill(5)+"-"+str(int(fiberid)).zfill(4)+".fits" )
get_path_to_spectrum_26 = lambda plate, mjd, fiberid : os.path.join(os.environ['HOME'], 'SDSS', '26', 'spectra', str(int(plate)).zfill(4), "spec-"+str(int(plate)).zfill(4)+"-"+str(int(mjd)).zfill(5)+"-"+str(int(fiberid)).zfill(4)+".fits" )

class SpectraStackingEBOSS:
	"""
	The model luminosity function class
	:param in_file: file containing spectra ids to be stacked
	:param Resolution: Resolution
	:param out_file: where to output stacks
	"""
	def __init__(self, in_file, out_file, dLambda = 0.0001, dV=-9999.99):
		print( "input list:", in_file )
		self.in_file = in_file
		self.plates, self.mjds, self.fiberids, self.redshifts = n.loadtxt(self.in_file, unpack=True)
		self.out_file = out_file
		self.dLambda = dLambda
		self.wave= 10**n.arange(2.6, 4.0211892990699383, dLambda) # 1500,10500
		self.R = int(1/n.mean((self.wave[1:] -self.wave[:-1])/ self.wave[1:]))
		print( "R=", n.median(self.R) )
		self.dV = dV
		self.survey="eBOSS"
		#self.run2d = run2d
		#self.run1d = self.run2d
		#self.topdirBOSS = os.path.join(os.environ['BOSS_SPECTRO_REDUX'], run2d)

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

	def getSpectra(self, path_to_spectrum):
		hdulist = fits.open(path_to_spectrum)
		wave = 10**hdulist[1].data['loglam']
		flux = hdulist[1].data['flux']
		ivar = hdulist[1].data['ivar']
		ratio = n.min(abs(10000.*n.log10(n.outer(wave, 1./maskLambda))), axis=1)
		margin = 1.5
		veto_sky = ratio <= margin
		selection = (veto_sky) & (ivar<=0) & (flux<0.)& (n.isinf(ivar)) & (n.isinf(flux))
		flux[selection] = n.zeros_like(ivar[selection])
		ivar[selection] = n.zeros_like(ivar[selection])
		out_sel = (flux>0)&(ivar>0)
		self.fluxl =flux[out_sel]
		self.fluxlErr=ivar[out_sel]**(-0.5)
		self.wavelength = wave[out_sel] #/(1+z)

	def stackSpectra(self):
		"""
		Function that constructs the stacks for a luminosity function. 
		It loops over the list of spectra given in the catalog of the LF. 
		First it sorts the catalog by the line luminosity. 
		And then stacks the first Nspec, then the next Nspec together.
		"""
		# loop over the file with N sorted with luminosity
		specMatrix, specMatrixErr, specMatrixWeight=[],[],[]
		
		for plate, mjd, fiber, redshift in zip(self.plates, self.mjds, self.fiberids, self.redshifts):
			try:
				#print(plate, mjd, fiber, redshift)
				if plate > 3006 :
					path_to_spectrum = get_path_to_spectrum_v5_10_0(plate, mjd, fiber)
				else:
					path_to_spectrum = get_path_to_spectrum_26(plate, mjd, fiber)
					
				if os.path.isfile(path_to_spectrum):
					self.getSpectra(path_to_spectrum)
					pts,ptsErr = self.convertSpectrum(redshift)
					specMatrix.append(pts)
					specMatrixErr.append(ptsErr)
					weight=1.
					specMatrixWeight.append(n.ones_like(pts)*weight)
				else: # for ELG spectra in v5_10_7
					path_to_spectrum = get_path_to_spectrum_v5_10_7(plate, mjd, fiber)
					if os.path.isfile(path_to_spectrum):
						self.getSpectra(path_to_spectrum)
						pts,ptsErr = self.convertSpectrum(redshift)
						specMatrix.append(pts)
						specMatrixErr.append(ptsErr)
						weight=1.
						specMatrixWeight.append(n.ones_like(pts)*weight)
			except(ValueError):
				print('value error !')

		specMatrixWeight=n.array(specMatrixWeight)
		specMatrix=n.array(specMatrix)
		specMatrixErr=n.array(specMatrixErr)
		print( "now stacks" )
		wavelength, medianStack, meanStack, meanWeightedStack, jackknifStackErrors, jackknifeSpectra, NspectraPerPixel = self.stack_function( specMatrix ,specMatrixWeight)
		cols = fits.ColDefs([wavelength, medianStack, meanStack, meanWeightedStack, jackknifStackErrors, jackknifeSpectra, NspectraPerPixel])
		tbhdu = fits.BinTableHDU.from_columns(cols)
		prihdr = fits.Header()
		prihdr['author'] = "JC"
		prihdr['survey'] = self.survey
		prihdr['in_file'] = os.path.basename(self.in_file)[:-4]
		prihdr['Nspec'] = len(self.plates)
		prihdu = fits.PrimaryHDU(header=prihdr)
		thdulist = fits.HDUList([prihdu, tbhdu])
		if os.path.isfile(self.out_file):
			os.remove(self.out_file)
		print( "stack written to", self.out_file )
		thdulist.writeto(self.out_file)

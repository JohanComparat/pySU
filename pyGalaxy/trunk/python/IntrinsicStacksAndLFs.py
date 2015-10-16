"""
*Imports*::

	import os 
    from os.path import join
	import astropy.cosmology as co
	cosmo=co.FlatLambdaCDM(H0=70,Om0=0.3)
	import astropy.io.fits as fits
	import numpy as n

"""

import os 
from os.path import join

import astropy.cosmology as co
cosmo=co.FlatLambdaCDM(H0=70,Om0=0.3)
import astropy.io.fits as fits
import numpy as n
import glob
from scipy.interpolate import interp1d


# /database/Products_Galaxies/emissionLineLuminosityFunctions/H1_4862/
# H1_4862_z0300.fitInfo
# H1_4862_z0300.points

class IntrinsicStacksAndLFs:

	"""
	This class reconnects the stacks and the luminosity function together
	:param LFfile: fits file of the catalog of the LF
	:param fireflyModel: firefly library used

	:param lineName: name of the line used in the catalogs.
	:param cosmology: cosmology used (astropy class) Default H0=70,Omega matter=0.3
    :param surveyName: Name of the survey used (needs to be the one given in the database)
    :param redshift_catalog: name of the redshift catalog 
	:param SNlimit: signal to noise ratio limit. For a line measurement to be included SN of the line must be greater than SNlimit.
	:param luminosityBins: bins in luminosity equally spaced in log space.
	:param Nstack: number of spectra per stack
	:param Nclustering:	number of spectra per clustering sample from the brightest to the faintest
	:param outputFolder: folder where the results will be written
	:param zmin: minimum redshift included
	:param zmax: maximum redshift included
	"""
	def __init__(self, LF_file="/home/comparat/database/DEEP2/products/emissionLineLuminosityFunctions/H1_4862/H1_4862-DEEP2-z0.746.fits", fireflyModel = "MarastonUVext"
):
		self.LF_file = LF_file
		self.stackList = n.array(glob.glob( join("/home/comparat/database/LFstacks/" , self.LF_file.split('/')[-1][:-5] + "_stack_*"+fireflyModel +"-modeled.fits") ))
		self.LF_measurement = self.LF_file[:-5] + ".txt"

		self.survey=LF_file.split('/')[-1].split('-')[1]
		self.line=LF_file.split('/')[-1].split('-')[0]
		self.lineDict = {'O2_3728' : r'$[O^{3728}_{II}]$', 'O3_5007' : r'$[O^{5007}_{III}]$', 'H1_4862' : r'$H^{4861}_{\beta}$', 'H1_6564': r'$H^{6564}_{\alpha}$'}
		self.lineLabel = self.lineDict[self.line]
		self.lineDictCorr4341 = {'O2_3728' : 'EBV_4862_4341_CORRO2', 'O3_5007' : 'EBV_4862_4341_CORRO3', 'H1_4862' : 'EBV_4862_4341_CORRHb'}
		self.lineDictCorr4341Err = {'O2_3728' : 'EBV_4862_4341_CORRO2_err', 'O3_5007' : 'EBV_4862_4341_CORRO3_err', 'H1_4862' : 'EBV_4862_4341_CORRHb_err'}
		self.lineCorr4341Name = self.lineDictCorr4341[self.line]
		self.lineCorr4341ErrName = self.lineDictCorr4341Err[self.line]

		self.lineDictCorr4102 = {'O2_3728' : 'EBV_4862_4102_CORRO2', 'O3_5007' : 'EBV_4862_4102_CORRO3', 'H1_4862' : 'EBV_4862_4102_CORRHb'}
		self.lineDictCorr4102Err = {'O2_3728' : 'EBV_4862_4102_CORRO2_err', 'O3_5007' : 'EBV_4862_4102_CORRO3_err', 'H1_4862' : 'EBV_4862_4102_CORRHb_err'}
		self.lineCorr4102Name = self.lineDictCorr4102[self.line]
		self.lineCorr4102ErrName = self.lineDictCorr4102Err[self.line]

		self.catalog=fits.open(self.LF_file)
		self.completeness = self.catalog[0].header['COMPLETENESS']
		self.volume = self.catalog[0].header['VOLUME']

		

	def define_correction(self):
		correction=[]
		for el in self.stackList:
			names = el.split('_')
			luminosity = float(names[-2])
			if luminosity > self.completeness:
				stackModel=fits.open(el)
				#print el, stackModel[0].header
				try:
					correction.append([luminosity, stackModel[0].header["flux_" + self.line + "_intrinsic"]/ stackModel[0].header[self.line + "_flux_nc"], stackModel[0].header["flux_" + self.line + "_intrinsic_err"]/ stackModel[0].header["flux_" + self.line + "_intrinsic"] ])		
				except KeyError :
					pass

		self.correction = n.transpose(correction)
		print self.correction

	def apply_correction(self):
		self.define_correction()
		if len(self.correction)>0:
			Lmin, Lmax, Lmean, phi, phiErr, phiErr_poisson, ngals = 	n.loadtxt( self.LF_measurement,unpack=True)
			x=n.hstack(( n.min(Lmin),self.correction[0],n.max(Lmax) ))
			y=n.hstack(( self.correction[1][0],self.correction[1], self.correction[1][-1] ))
			yErr=n.hstack(( self.correction[2][0],self.correction[2], self.correction[2][-1] ))
			corr = interp1d(x,y)	

			Lmin_c = Lmin[1:-1] / corr(Lmin[1:-1]) 
			Lmax_c = Lmax[1:-1] / corr(Lmax[1:-1]) 
			Lmean_c = Lmean[1:-1] / corr(Lmean[1:-1]) 

			f=open(self.LF_measurement[:-4] + "_intrinsic.txt",'w')
			n.savetxt(f, n.transpose([Lmin_c, Lmax_c, Lmean_c, phi[1:-1], phiErr[1:-1], phiErr_poisson[1:-1], ngals[1:-1]]), header = " Lmin Lmax Lmean phi phiErr phiErr_poisson Ngalaxy")
			f.close()

			f=open(self.LF_measurement[:-4] + "_correction.txt",'w')
			n.savetxt(f, n.transpose(self.correction), header = " L EBV_4862_4341_CORRHb EBV_4862_4341_CORRHb_err EBV_4862_4102_CORRHb EBV_4862_4102_CORRHb ")
			f.close()


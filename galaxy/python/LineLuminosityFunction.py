"""
.. class:: LineLuminosityFunction
.. class:: ModelLuminosityFunction

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class LineLuminosityFunction is dedicated to measuring the line luminosity functions. The class ModelLuminosityFunction is dedicated to fitting models to the LFs.

"""
import os 
from os.path import join
import cPickle
import astropy.cosmology as co
cosmo=co.Planck15 # co.FlatLambdaCDM(H0=70,Om0=0.3)
import astropy.io.fits as fits
import numpy as n
from scipy.optimize import curve_fit

class ModelLuminosityFunction:
	"""
	The model luminosity function class
	:param lineWavelength: restframe wavelength in the air
	:param lineName: name of the line used in the catalogs.
	:param cosmology: cosmology used (astropy class) Default H0=70,Omega matter=0.3
	:param LF_file_list: list of LF_files. Each entry is a list of 2 files, the fits file (data used to derivee the LF) and the txt file (LF).
	:param model: Model to be fitted : Saunders, Saunders3P, or Schechter, DoublePL
	:param p0: Initial guess of the parameters
	:param outputDirectory: where to output the fits
	:param fileName: file name where things are saved.
	"""
	def __init__(self, lineWavelength=3727.4228417998916, lineName="O2_3728", cosmology = cosmo, LF_file_list=[], model='Saunders', p0=[1.], outputDirectory="", fileName="", fixedSigma = 0.5, fixedAlpha = -1.5):
		self.lineWavelength = lineWavelength
		self.lineName = lineName
		self.cosmology = cosmology
		self.database_dir = os.environ['DATA_DIR']
		self.biblioPts_dir = join(self.database_dir , "Products_Galaxies","emissionLineLuminosityFunctions","bibliographyPoints")
		self.fixedSigma = fixedSigma
		self.fixedAlpha = fixedAlpha

		self.LF_file_list = LF_file_list

		self.schecterFct=lambda logl,logls,ps,a : ps * (10**logl/10**logls)**(a+1) * n.e**(-10**logl/10**logls)
		self.saundersFct=lambda logl,logls,ps,a,sig : ps * (10**logl/10**logls)**(a+1) * n.e**( -n.log10( 1 +10**logl/10**logls)**2./(2*sig**2.))
		self.saundersFct3P=lambda logl,logls,ps,a : self.saundersFct(logl, logls, ps, a, self.fixedSigma)
		self.saundersFct2P=lambda logl,logls,ps : self.saundersFct(logl, logls, ps, self.fixedAlpha, self.fixedSigma)
		self.doublePL=lambda logl,a,b,logls : 10**( (1+a)*(logl - logls) + b )
		modelDict = {'Saunders2P':self.saundersFct2P, 'Saunders3P':self.saundersFct3P, 'Saunders': self.saundersFct, 'Schechter': self.schecterFct, 'DoublePL' : self.doublePL }
		self.model = model
		self.model_to_fit = modelDict[self.model]
		self.p0 = p0
		self.outputDirectory = outputDirectory
		self.fileName = fileName

	def fitModel(self,bibliographyPts=0):
		"""
		fits the chosen model to the list of Measured LFs
		"""
		x,y,ye,name=[],[],[],[]
		for el in self.LF_file_list:
			lf_fits_file,lf_measurement_file = el

			hduG=fits.open(lf_fits_file)
			hdu=hduG[1]
			catalog=hdu.data
			completeness=hduG[0].header['COMPLETENESS']
			volume=hduG[0].header['VOLUME']
			Lmin, Lmax, Lmean, phi, phiErr, phiErr_poisson, ngals= n.loadtxt( lf_measurement_file, unpack=True)
			selected = n.array([LmeanEl > completeness for LmeanEl in Lmean]) & (phi > 0) & (phi > 5. / volume)
			x.append(Lmean[selected])
			y.append(phi[selected])
			ye.append(phiErr[selected])
			name.append(lf_measurement_file.split('/')[-1].split('-')[1])
			print lf_fits_file, completeness, Lmean[selected], phi[selected]

		if bibliographyPts != 0 :
			for jj in range(len(bibliographyPts)):
				x.append(bibliographyPts[jj][0])
				y.append(bibliographyPts[jj][1])
				ye.append(bibliographyPts[jj][2])
				name.append(bibliographyPts[jj][3])
			
		xFi=n.hstack((x))
		yFi=n.hstack((y))
		yeFi=n.hstack((ye))

		ok=(xFi>0)&(yFi>0)&(yeFi>0)
		# print xFi,yFi,yeFi

		xF=xFi[ok][n.argsort(xFi[ok])]
		yF=yFi[ok][n.argsort(xFi[ok])]
		yeF=yeFi[ok][n.argsort(xFi[ok])]
		print "final array", len(xF), len(yF), len(yeF)
		print xF, yF, yeF
		popt,popc=curve_fit(self.model_to_fit,n.log10(xF),yF,p0=self.p0, sigma=yeF, maxfev=5000000)
		print "result", popt, popc
		#yMod=self.model_to_fit(xF,popt)
		n.savetxt(join(self.outputDirectory,self.fileName+".points"),n.transpose([xF,yF,yeF]))
		f = open(join(self.outputDirectory,self.fileName+".fitInfo"),'w')
		for el in self.LF_file_list :
			f.write(el[0])
			f.write(' \n')
			f.write(el[1])
			f.write(' \n')

		f.write(self.model)
		f.write(' \n')
		f.write(str(popt))
		f.write(' \n')
		f.write(str(popc))
		f.write(' \n')
		f.close()
		if len(popt)==4:
			f = open(join(self.outputDirectory,self.fileName+".fitInfo.PKL"),'w')
			cPickle.dump([popt,popc,self.fileName],f)
			f.close()
		if len(popt)==3:
			f = open(join(self.outputDirectory,self.fileName+".fitInfo.PKL"),'w')
			cPickle.dump([popt,popc,self.fileName, self.fixedSigma],f)
			f.close()
		if len(popt)==2:
			f = open(join(self.outputDirectory,self.fileName+".fitInfo.PKL"),'w')
			cPickle.dump([popt,popc,self.fileName, self.fixedSigma, self.fixedAlpha],f)
			f.close()
		return popt,popc,xF,yF,yeF, x, y, ye, name



class LineLuminosityFunction:
	"""
	The line luminosity function class
	:param lineWavelength: restframe wavelength in the air
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
	def __init__(self, lineWavelength=3727.4228417998916, lineName="O2_3728", cosmology = cosmo, surveyName =  "DEEP2", redshift_catalog = "zcat.deep2.dr4.v2.LFcatalog.fits", SNlimit = 5., luminosityBins = n.logspace(38,45,50), Nstack = 400, Nclustering = 400, outputFolder="emissionLineLuminosityFunctions/" , zmin=0.6, zmax=0.8):
		self.lineWavelength = lineWavelength
		self.lineName = lineName
		self.cosmology = cosmology
		self.surveyName = surveyName
		self.redshift_catalog = redshift_catalog
		self.database_dir = os.environ['DATA_DIR']
		self.survey_dir = join( self.database_dir, self.surveyName)
		self.catalog_dir = join(self.survey_dir,"catalogs")
		self.output_dir = join(self.survey_dir,"products",outputFolder,lineName)
		os.system('mkdir '+self.output_dir)
		hd = fits.open(join(self.catalog_dir,self.redshift_catalog))
		self.catalog = hd[1].data
		hd.close()
		self.Ngalaxies = len(self.catalog)
		self.SNlimit =  SNlimit # 5.
		#self.nbins = 15#n.arange(38.5,45,0.25)#15
		self.luminosityBins = luminosityBins #15
		#self.nbinsUD = 4
		self.Nstack = Nstack # number of spectra per stack
		self.Nclustering = Nclustering # number of spectra per clustering sample from the brightest to the 

		self.zmin = zmin
		self.zmax = zmax

		self.luminosity = self.catalog[lineName+'_luminosity']
		self.luminosityErr = self.catalog[lineName+'_luminosityErr']
		self.volume_per_sq_degree=lambda z1,z2 : (cosmo.comoving_volume( z2 ) - cosmo.comoving_volume( z1 )) *n.pi/129600.

		self.lineSelection=(self.catalog[lineName+'_flux']>0.)& (self.catalog[lineName+'_fluxErr'] >0.) & (self.catalog[lineName+'_flux'] > self.SNlimit * self.catalog[lineName+'_fluxErr']) & (self.catalog[lineName+'_luminosity']>0)& (self.catalog[lineName+'_luminosity']<1e50)

	def setRedshiftArray(self,redshiftColumn='ZBEST'):
		""" sets the redshift array
		:param redshiftColumn: column of the catalog corresponding to the redshift.
		Stores it in self.redshift.
		"""
		self.redshift = self.catalog[redshiftColumn]

	def setWeightArray(self,weightColumn, area):
		""" sets the weight column 
		:param weightColumn: statistical weight per galaxy 1 / (area * TSR * SSR)
		Divides the weight by the volume of the bin stores it in self.weight.
		"""
		self.volume =  self.volume_per_sq_degree(self.zmin,self.zmax) * area
		self.weight = weightColumn / self.volume_per_sq_degree(self.zmin,self.zmax)

	def setRedshiftSelection(self, redshiftQualityColumn='ZQUALITY', upperBound=0.9, lowerBound=7.):
		""" sets the redshift selection
		:param redshiftQualityColumn: column of the catalog corresponding to the  quality of the redshifts.
		:param lowerBound : lower bound to redshift quality :  zquality > lowerBound
		:param upperBound : upper bound to the redshift quality : zquality < upperBound
		Stores it in self.redshiftSelection.
		"""
		self.redshiftSelection = ( self.redshift>self.zmin ) & ( self.redshift<self.zmax ) & (self.catalog[redshiftQualityColumn]> lowerBound) & (self.catalog[redshiftQualityColumn] < upperBound)

	def computeMeanWeightedRedshift(self,sel):
		""" Computes the weighted mean redshift of the sample.
		"""
		selection = (sel) & (self.lineSelection) & (self.redshiftSelection)
		self.meanRedshift = n.average(self.redshift[selection], weights = self.weight[selection])

	def computeHistogramLF(self,sel):
		""" Computes the weighted and unweighted histogram to get the number density and Poisson errors.
		:param sel: array selecting the galaxies of interest in the catalog (Boolean). 
		Returns Weighted density, Error on the weighted density, Number of galaxies used in eacah bin, the luminosity bins.
		It stores the values in self.LF, self.LFerr_poisson, self.ngals. It also evaluates the mean luminosity in each luminosity bin self.xL and dlogL to obtain the LF
		"""
		selection = (sel) & (self.lineSelection) & (self.redshiftSelection)
		N10p,bin1p=n.histogram(self.luminosity[selection],bins=self.luminosityBins)
		N10,bin1=n.histogram(self.luminosity[selection], bins= self.luminosityBins, weights= self.weight[selection] )
		self.LF, self.LFerr_poisson, self.ngals = N10, N10*N10p**0.5/N10p, N10p

		xSelections=n.array([ (self.luminosity > self.luminosityBins[ii]) &(self.luminosity< self.luminosityBins[ii+1] ) & (selection) for ii in range( len( self.luminosityBins ) -1 ) ])

		xLi= []
		for jj in range(len(xSelections)) :
			if len(self.luminosity[xSelections[jj]])>0:
				xLi.append( n.average( self.luminosity[xSelections[jj]], weights= self.weight[xSelections[jj]] ) )
			else:
				xLi.append( (self.luminosityBins[jj]+self.luminosityBins[jj+1])/2. )
		
		self.xL=n.array(xLi)
		dLogL_all = (self.luminosityBins[1:] - self.luminosityBins[:-1]) / ((self.luminosityBins[1:] + self.luminosityBins[:-1])/2.)
		self.dLogL = dLogL_all[0]

	def computeHistogramVariance(self,sel,jk=0.1):
		""" Computes the variance of the histogram using N subsamples.
		:param sel: array selecting the galaxies of interest in the catalog (Boolean). 
		:param jk: percentage of the data set removed in each realization.
		Stores the values in self.LFerr_jackknife
		"""
		selection = (sel) & (self.lineSelection) & (self.redshiftSelection)
		#N10p,bin1p=n.histogram(self.luminosity[selection],bins=self.luminosityBins)
		L_jk = self.luminosity[selection]
		w_jk = self.weight[selection]
		rdArr=n.random.rand(len(L_jk))
		values=n.arange(0,1+0.9*jk,jk)
		randSelNot=n.array([(rdArr>values[jj])&(rdArr<values[jj+1]) for jj in range(len(values)-1)])
		randSel=n.array([(el==False) for el in randSelNot])
		lumJK=[]
		for selR in randSel :
			N10,bin1=n.histogram(L_jk[selR], bins= self.luminosityBins, weights= w_jk[selR] )
			lumJK.append(N10)

		self.LFerr_jackknife = n.std(lumJK,axis=0)

	def get_completness_limit(self,sel):
		"""
		Estimates the completeness. It maps the maximum of the EW distribution to a Luminosity limit.
		"""
		selection = (sel) & (self.lineSelection) & (self.redshiftSelection) & (self.catalog[self.lineName+'_EW']!=n.inf)&(self.catalog[self.lineName+'_EW']>0)&(self.catalog[self.lineName+'_EW']<10**4)
		bins=n.logspace(0,3,40)
		print len(self.catalog[self.lineName+'_EW'][selection]), n.min(self.catalog[self.lineName+'_EW'][selection]), n.max(self.catalog[self.lineName+'_EW'][selection])
		aa,bb = n.histogram(self.catalog[self.lineName+'_EW'][selection], bins=bins)
		self.completness_limit_EW = bb[n.argmax(aa)+1]
		
		EWselection = (self.catalog[self.lineName+'_EW'][selection] >0.9* self.completness_limit_EW )&( self.catalog[self.lineName+'_EW'][selection]<1.1* self.completness_limit_EW)

		self.completness_limit_luminosity = n.median( self.catalog[ self.lineName+'_luminosity'][ selection ][ EWselection ])

		# bins=n.logspace(39.5,43,20)
		# aa,bb = n.histogram(self.catalog[self.lineName+'_luminosity'][selection], bins=bins)
		#self.completness_limit_luminosity = bb[n.argmax(aa)+1]

	def correctVolumeEffect(self, mCorr, error):
		""" corrects from the volume effect on the measurement and the error estimation
		"""
		self.phi = self.LF/self.dLogL * mCorr
		self.phiErr = ( ( self.LFerr_jackknife / self.dLogL * mCorr)**2. + (self.phi * error )**2. )**0.5
		self.phiErr_poisson = self.LFerr_poisson/self.dLogL * mCorr

	def writeLF(self,sel,surveyNameSuffix=""):
		""" writes the measured LF and the data used to derive it to an ascii and a fits file.
		"""
		filename = self.lineName + "-" + self.surveyName+surveyNameSuffix + "-z" + str( n.round( self.meanRedshift ,3 )) 

		selection = (sel) & (self.lineSelection) & (self.redshiftSelection)
		new_columns = self.catalog.columns

		tbhdu = fits.BinTableHDU.from_columns(new_columns)
		tbhdu.data = tbhdu.data[selection]
		prihdr = fits.Header()
		prihdr['HIERARCH COMPLETENESS'] = self.completness_limit_luminosity
		print "volume",self.volume.value
		prihdr['VOLUME'] = self.volume.value
		prihdu = fits.PrimaryHDU(header=prihdr)
		thdulist = fits.HDUList([prihdu, tbhdu])

		os.system('rm -rf '+join(self.output_dir , filename + ".fits"))
		thdulist.writeto(join(self.output_dir , filename + ".fits"))                        

		head= "z_mean z_min z_max L_mean L_min L_max phi_mean phi_min phi_max  Ngalaxy"
		f=open(join(self.output_dir , filename + ".txt"),'w')
		n.savetxt(f, n.transpose([self.meanRedshift*n.ones_like(self.phi) , self.zmin*n.ones_like(self.phi) , self.zmax*n.ones_like(self.phi), self.xL, self.luminosityBins[:-1], self.luminosityBins[1:], self.phi, self.phi-self.phiErr, self.phi + self.phiErr, self.ngals]), header= head)
		f.close()

		head= " Lmin & Lmax  & Lmean & phi & phiErr & Ngalaxy"
		f=open(join(self.output_dir , filename + ".tex"),'w')
		n.savetxt(f, n.transpose([n.log10(self.luminosityBins[:-1][self.xL>self.completness_limit_luminosity]), n.log10(self.luminosityBins[1:][self.xL>self.completness_limit_luminosity]), n.log10(self.xL[self.xL>self.completness_limit_luminosity]), n.log10(self.phi[self.xL>self.completness_limit_luminosity]), n.log10(self.phiErr[self.xL>self.completness_limit_luminosity]), self.ngals[self.xL>self.completness_limit_luminosity]]) ,header= head, delimiter = "&", fmt='%10.2f', newline= " \\\\ ")
		f.close()

		head= " zMean & zmin & zmax  & logV & area & NgalaxyTotal & Lmin & Lcompleteness & Lmax & "+self.surveyName
		f=open(join(self.output_dir , filename + "summaryLine.tex"),'w')
		n.savetxt(f,  n.transpose([[self.meanRedshift] , [self.zmin] , [self.zmax] , [n.log10(self.volume.value)] , [self.area] , [n.sum(self.ngals[self.xL>self.completness_limit_luminosity]).T] , [n.log10(n.min(self.luminosity[selection]))] , [n.log10(self.completness_limit_luminosity)], [n.log10(n.max(self.luminosity[selection]))] ]) ,header= head, delimiter = "&", fmt='%10.2f', newline= " \\\\ ")
		f.close()


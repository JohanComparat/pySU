"""
.. class:: LineFittingLibrary

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

This class contains a variety of function to fit emission or absorption lines in galaxy spectra.

"""
from scipy.optimize import curve_fit
import numpy as n
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p
from scipy.interpolate import interp1d
from scipy.integrate import quad
# Location of the emission lines of interest:
import astropy.constants as cc
c=cc.c.value # speed of light
#from lineList import *

class LineFittingLibrary:
	"""
	Loads the environement proper to fit lines :
	 * Gaussian line model
	 * Lorentzian line model
	 * pseudoVoigt line model
	 * conversion magnitude AB to flux : flambda to fnu

	:param dV: the default value  (def: -9999.99)
	:param fitWidth: width in Angstrom around the line where the fit is performed, default 35 Angstrom 
	"""
	def __init__(self,dV=-9999.99):
		self.dV=dV # default value put in the catalogs
		# Line models
		self.gaussianLine=lambda aa,sigma,F0,a0,continu : continu + F0*(n.e**( -(aa-a0)**2. / (2.*sigma**2.)))/ (abs(sigma)*(2.*n.pi)**0.5)
		self.gaussianLineNC=lambda aa,sigma,F0,a0 : F0*(n.e**(-(aa-a0)**2./ (2.*sigma**2.) ))/(abs(sigma)*(2.*n.pi)**0.5)
		self.lorentzLine=lambda aa,gamma,F0,a0,continu  : continu + F0 * abs(gamma) / (n.pi* ((aa-a0)**2 +gamma**2))
		self.pseudoVoigtLine=lambda aa,fwhm,F0,a0,continu,sh : continu + F0*abs(sh)/(1+ ((aa-a0) /(fwhm/2.))**2.)+F0*(1-abs(sh))*n.e**( -n.log(2)* ((aa-a0)/(fwhm/2.))**2.) 

		# conversion magnitude flux
		self.fnu = lambda mAB : 10**(-(mAB+48.6)/2.5) # erg/cm2/s/Hz
		self.flambda= lambda mAB, ll : 10**10 * c * self.fnu(mAB) / ll**2. # erg/cm2/s/A

	def integrateMAG(self,wl,spec1d,err1d,filt,xmin=5000.,xmax=7500.):
		"""
		Integrates a spectrum over a filter curve.

		:param wl: wavelength (array)
		:param spec1d: flux, f lambda convention (array)
		:param err1d: flux error (array)
		:param filt: filter curve (interpolation 1d)
		:param xmin: lower integration boundary (Angstrom)
		:param xmax: higher integration boundary (Angstrom)

		returns :
		 * integral of filter curve
		 * integral of spec1d
		 * integral of spec1d * filter curve
		 * integral of (spec1d + err1d) * filter curve
		 * integral of (spec1d - err1d) * filter curve
	 
		"""
		filtTp=filt(wl)
		Lfilt=quad(filt,xmin,xmax,limit=500000)[0]
		toInt=interp1d(wl,spec1d)
		Lspec=quad(toInt,xmin,xmax,limit=500000)[0]
		toInt=interp1d(wl,spec1d*filtTp)
		Lg=quad(toInt,xmin,xmax,limit=500000)[0]
		toInt=interp1d(wl,(spec1d+err1d)*filtTp)
		LgU=quad(toInt,xmin,xmax,limit=500000)[0]
		toInt=interp1d(wl,(spec1d-err1d)*filtTp)
		LgL=quad(toInt,xmin,xmax,limit=500000)[0] 
		return Lfilt, Lspec, Lg, LgU, LgL

	def getFractionObsMed(self,mag,lambdaMag,fl,flErr):
		"""
		Computes the fraction of light captured by the spectrograph in a broad band by comparing the median flux in the broad band to the magnitude converted to flux at the mean wavelength of the broad band.

		:param mag: magnitude AB (float, mag)
		:param lambdaMag: mean wavelength covered by the magnitude AB (float, Angstrom)
		:param fl: flux observed in the broad band (array, f lambda)
		:param flErr: error on the flux observed in the broad band (array, f lambda)

		Returns	:
		 * fraction of light observed
		 * error on the fraction of light observed
		"""
		goal=self.flambda(mag,lambdaMag)
		fo=goal/n.median(fl)
		fo_err=goal/n.median(flErr)
		return fo, fo_err

	def getFractionObsMag(self,mag,lambdaMag,filter,xmin,xmax,wl,fl,flErr):
		"""
		Computes the fraction of light captured by the spectrograph in a broad band by comparing the integrated flux in the broad band to the magnitude.

		:param mag: magnitude AB (float, mag)
		:param lambdaMag: mean wavelength covered by the magnitude AB (float, Angstrom)
		:param fl: flux observed in the broad band (array, f lambda)
		:param flErr: error on the flux observed in the broad band (array, f lambda)
		:param filt: filter curve (interpolation 1d)
		:param xmin: lower integration boundary (Angstrom)
		:param xmax: higher integration boundary (Angstrom)

		Returns	:
		 * fraction of light observed
		 * error on the fraction of light observed
		"""
		goal=self.flambda(mag,lambdaMag)
		Lfilt, Lspec, Lg, LgU, LgL=self.integrateMAG(wl,fl,flErr,filter,xmin,xmax)
		fo=Lg/Lfilt/goal
		fo_err=(LgU/Lfilt/goal-LgL/Lfilt/goal)/2
		return fo, fo_err

	def plotLineFit(self,wl,fl,flErr,lineModel,a0,path_to_fig="plot.pdf", title=" - "):
		"""
		Plots a spectrum and the emission line model fitted.

		:param wl: wavelength (array, Angstrom)
		:param fl: flux observed in the broad band (array, f lambda)
		:param flErr: error on the flux observed in the broad band (array, f lambda)
		:param lineModel: model output by the line fitting functions (array, (2,N) wavelength and flux)
		:param a0: position of the peak of the line
		:param path_to_fig: where you wish to save the figure
		"""
		p.figure(0,(8,4))
		p.plot(wl,fl,'k')
		p.plot(wl,fl+flErr,'g--')
		p.plot(wl,fl-flErr,'g--')
		p.axvline(a0)
		wlMod=lineModel[0]
		p.plot(wlMod,lineModel[1],'r')
		p.xlim((wlMod.min()-50,wlMod.max()+50))
		p.yscale('log')
		p.ylim((n.max([lineModel[1].min() / 5., 1e-18]), lineModel[1].max() * 5.))
		p.title(title)
		p.savefig(path_to_fig)
		p.clf()

	def fit_Line_position_C0noise(self,wl,spec1d,err1d,a0=5007.,lineName="AL",fitWidth=20,DLC=20, p0_sigma=15.,p0_flux=8e-17,p0_share=0.5,continuumSide="left",model="gaussian"):
		"""
		fits a line profile to a spectrum around a fixed line position

		:param wl: wavelength (array, Angstrom)
		:param spec1d: flux observed in the broad band (array, f lambda)
		:param err1d: error on the flux observed in the broad band (array, f lambda)
		:param a0: expected position of the peak of the line in the observed frame (redshifted). a0 is not fitted, it is given.
		:param lineName: suffix characterizing the line in the headers of the output
		:param DLC: wavelength extent to fit the continuum around the line. (def: 230 Angstrom)
		:param p0_sigma: prior on the line width in A (def: 15 A)
		:param p0_flux: prior on the line flux in erg/cm2/s/A (def: 8e-17)
		:param p0_share: prior on the share of Gaussian and Lorentzian model. Only used if the line is fitted with a pseudoVoigt profile width (def: 0.5 no units)
		:param continuumSide: "left" = bluewards of the line or "right" = redwards of the line
		:param model: line model to be fitted : "gaussian", "lorentz" or "pseudoVoigt".

		Returns :
		 * array 1 with the parameters of the model
		 * array 2 with the model (wavelength, flux model)
		 * header corresponding to the array 1
		"""
		header=" "+lineName+"_a0 "+lineName+"_flux "+lineName+"_fluxErr "+lineName+"_sigma "+lineName+"_sigmaErr "+lineName+"_continu "+lineName+"_continuErr "+lineName+"_EW "+lineName+"_fd_a0_l "+lineName+"_fd_a0_r "+lineName+"_chi2 "+lineName+"_ndof"
		headerPV=" "+lineName+"_a0 "+lineName+"_flux "+lineName+"_fluxErr "+lineName+"_sigma "+lineName+"_sigmaErr "+lineName+"_continu "+lineName+"_continuErr "+lineName+"_EW "+lineName+"_share "+lineName+"_shareErr "+lineName+"_fd_a0_l "+lineName+"_fd_a0_r "+lineName+"_chi2 "+lineName+"_ndof"
		outPutNF=n.array([a0, self.dV,self.dV,self.dV, self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV])
		outPutNF_PV=n.array([a0, self.dV,self.dV,self.dV, self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV])
		modNF=n.array([self.dV,self.dV])
		if continuumSide=="left":
			domainLine=(wl>a0-fitWidth)&(wl<a0+fitWidth)
			domainCont=(wl>a0-DLC-fitWidth)&(wl<a0-fitWidth)
			if a0<wl.max()-DLC and a0>wl.min()+fitWidth and a0<wl.max()-fitWidth and len(domainLine.nonzero()[0])>2 and len(domainCont.nonzero()[0])>2 :
				continu=n.median(spec1d[domainCont])
				continuErr=n.median(err1d[domainCont])
				if model=="gaussian":
					flMod=lambda aa,sigma,F0,a0,continu : self.gaussianLine(aa,sigma,F0,a0,continu)
					p0=n.array([p0_sigma,p0_flux,a0,continu])
				if model=="lorentz":
					flMod=lambda aa,sigma,F0,a0,continu : self.lorentzLine(aa,sigma,F0,a0,continu)
					p0=n.array([p0_sigma,p0_flux, a0,continu])
				if model=="pseudoVoigt":
					flMod=lambda aa,sigma,F0,sh,a0,continu : self.pseudoVoigtLine(aa,sigma,F0,a0,continu,sh)
					p0=n.array([p0_sigma,p0_flux,p0_share,a0,continu])
				interp=interp1d(wl,spec1d)
				fd_a0_r=interp(a0+0.2)
				fd_a0_l=interp(a0-0.2)
				if fd_a0_r>continu and fd_a0_l>continu :
					out = curve_fit(flMod, wl[domainLine], spec1d[domainLine], p0=p0,sigma=continu*n.ones_like(err1d[domainLine]),maxfev=1000000000, gtol=1.49012e-8)
					if out[1].__class__==n.ndarray and ( model=="gaussian" or model=="lorentz") : 
						model1=flMod(wl[domainLine],out[0][0],out[0][1],out[0][2],out[0][3])
						var=continu*n.ones_like(err1d[domainLine])
						chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
						ndof=len(var)
						sigma=out[0][0]
						sigmaErr=out[1][0][0]**0.5
						flux=out[0][1]
						fluxErr=out[1][1][1]**0.5
						a0=out[0][2]
						a0_err=out[1][2][2]**0.5
						continu=out[0][3]
						continuErr=out[1][3][3]**0.5
						EW=flux/continu
						outPut=n.array([a0,flux,fluxErr,sigma,sigmaErr,continu,continuErr,EW,fd_a0_l,fd_a0_r,chi2,ndof ])
						mod=n.array([wl[domainLine],model1])
						return outPut,mod,header
					elif model=="gaussian" or model=="lorentz" :
						return n.array([a0,self.dV,self.dV,self.dV, self.dV,continu,continuErr,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV ]),modNF,header
					elif out[1].__class__==n.ndarray and model=="pseudoVoigt" : 
						model1=flMod(wl[domainLine],out[0][0],out[0][1],out[0][2],out[0][3],out[0][4])
						var=continu*n.ones_like(err1d[domainLine])
						chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
						ndof=len(var)
						sigma=out[0][0]
						sigmaErr=out[1][0][0]**0.5
						flux=out[0][1]
						fluxErr=out[1][1][1]**0.5
						share=out[0][2]
						shareErr=out[1][2][2]**0.5
						a0=out[0][3]
						a0_err=out[1][3][3]**0.5
						continu=out[0][4]
						continuErr=out[1][4][4]**0.5
						EW=flux/continu
						outPut=n.array([a0,flux,fluxErr,sigma,sigmaErr,continu,continuErr,EW,share,shareErr,fd_a0_l,fd_a0_r,chi2,ndof])
						mod=n.array([wl[domainLine],model1])
						return outPut,mod,headerPV
					else :
						return n.array([a0,self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,self.dV,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,headerPV
				else :
					if  model=="gaussian" or model=="lorentz" :
						return n.array([a0,self.dV,self.dV,self.dV, self.dV,continu,continuErr,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV ]),modNF,header
					if model=="pseudoVoigt" :
						return n.array([a0,self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,self.dV,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,headerPV
			else :
				print "not enough space to fit the line"
				if  model=="gaussian" or model=="lorentz" :
					return outPutNF,modNF,header
				if model=="pseudoVoigt" :
					return outPutNF_PV,modNF,headerPV

		elif continuumSide=="right" :
			domainLine=(wl>a0-fitWidth)&(wl<a0+fitWidth)
			domainCont=(wl>a0+fitWidth)&(wl<a0+DLC+fitWidth)
			if a0<wl.max()-DLC and a0>wl.min()+fitWidth and a0<wl.max()-fitWidth and len(domainLine.nonzero()[0])>2 and len(domainCont.nonzero()[0])>2 :
				continu=n.median(spec1d[domainCont])
				continuErr=n.median(err1d[domainCont])
				if model=="gaussian":
					flMod=lambda aa,sigma,F0,a0,continu : self.gaussianLine(aa,sigma,F0,a0,continu)
					p0=n.array([p0_sigma,p0_flux,a0,continu])
				if model=="lorentz":
					flMod=lambda aa,sigma,F0,a0,continu : self.lorentzLine(aa,sigma,F0,a0,continu)
					p0=n.array([p0_sigma,p0_flux, a0,continu])
				if model=="pseudoVoigt":
					flMod=lambda aa,sigma,F0,sh,a0,continu : self.pseudoVoigtLine(aa,sigma,F0,a0,continu,sh)
					p0=n.array([p0_sigma,p0_flux,p0_share,a0,continu])
				interp=interp1d(wl,spec1d)
				fd_a0_r=interp(a0+0.2)
				fd_a0_l=interp(a0-0.2)
				if fd_a0_r>continu and fd_a0_l>continu :
					out = curve_fit(flMod, wl[domainLine], spec1d[domainLine], p0=p0,sigma=continu*n.ones_like(err1d[domainLine]),maxfev=1000000000, gtol=1.49012e-8)
					if out[1].__class__==n.ndarray and ( model=="gaussian" or model=="lorentz") : 
						model1=flMod(wl[domainLine],out[0][0],out[0][1],out[0][2],out[0][3])
						var=continu*n.ones_like(err1d[domainLine])
						chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
						ndof=len(var)
						sigma=out[0][0]
						sigmaErr=out[1][0][0]**0.5
						flux=out[0][1]
						fluxErr=out[1][1][1]**0.5
						a0=out[0][2]
						a0_err=out[1][2][2]**0.5
						continu=out[0][3]
						continuErr=out[1][3][3]**0.5
						EW=flux/continu
						outPut=n.array([a0,flux,fluxErr,sigma,sigmaErr,continu,continuErr,EW,fd_a0_l,fd_a0_r,chi2,ndof ])
						mod=n.array([wl[domainLine],model1])
						return outPut,mod,header
					elif model=="gaussian" or model=="lorentz" :
						return n.array([a0,self.dV,self.dV,self.dV, self.dV,continu,continuErr,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV ]),modNF,header
					elif out[1].__class__==n.ndarray and model=="pseudoVoigt" : 
						model1=flMod(wl[domainLine],out[0][0],out[0][1],out[0][2],out[0][3],out[0][4])
						var=continu*n.ones_like(err1d[domainLine])
						chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
						ndof=len(var)
						sigma=out[0][0]
						sigmaErr=out[1][0][0]**0.5
						flux=out[0][1]
						fluxErr=out[1][1][1]**0.5
						share=out[0][2]
						shareErr=out[1][2][2]**0.5
						a0=out[0][3]
						a0_err=out[1][3][3]**0.5
						continu=out[0][4]
						continuErr=out[1][4][4]**0.5
						EW=flux/continu
						outPut=n.array([a0,flux,fluxErr,sigma,sigmaErr,continu,continuErr,EW,share,shareErr,fd_a0_l,fd_a0_r,chi2,ndof])
						mod=n.array([wl[domainLine],model1])
						return outPut,mod,headerPV
					else :
						return n.array([a0,self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,self.dV,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,headerPV
				else :
					if  model=="gaussian" or model=="lorentz" :
						return n.array([a0,self.dV,self.dV,self.dV, self.dV,continu,continuErr,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV ]),modNF,header
					if model=="pseudoVoigt" :
						return n.array([a0,self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,self.dV,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,headerPV
			else :
				print "not enough space to fit the line"
				if  model=="gaussian" or model=="lorentz" :
					return outPutNF,modNF,header
				if model=="pseudoVoigt" :
					return outPutNF_PV,modNF,headerPV
					
	def fit_Line_position(self,wl,spec1d,err1d,a0=5007.,lineName="AL",fitWidth=20,DLC=20, p0_sigma=15.,p0_flux=8e-17,p0_share=0.5,continuumSide="left",model="gaussian"):
		"""
		fits a line profile to a spectrum around a fixed line position

		:param wl: wavelength (array, Angstrom)
		:param spec1d: flux observed in the broad band (array, f lambda)
		:param err1d: error on the flux observed in the broad band (array, f lambda)
		:param a0: expected position of the peak of the line in the observed frame (redshifted). a0 is not fitted, it is given.
		:param lineName: suffix characterizing the line in the headers of the output
		:param DLC: wavelength extent to fit the continuum around the line. (def: 230 Angstrom)
		:param p0_sigma: prior on the line width in A (def: 15 A)
		:param p0_flux: prior on the line flux in erg/cm2/s/A (def: 8e-17)
		:param p0_share: prior on the share of Gaussian and Lorentzian model. Only used if the line is fitted with a pseudoVoigt profile width (def: 0.5 no units)
		:param continuumSide: "left" = bluewards of the line or "right" = redwards of the line
		:param model: line model to be fitted : "gaussian", "lorentz" or "pseudoVoigt".

		Returns :
		 * array 1 with the parameters of the model
		 * array 2 with the model (wavelength, flux model)
		 * header corresponding to the array 1
		"""
		header=" "+lineName+"_a0 "+lineName+"_flux "+lineName+"_fluxErr "+lineName+"_sigma "+lineName+"_sigmaErr "+lineName+"_continu "+lineName+"_continuErr "+lineName+"_EW "+lineName+"_fd_a0_l "+lineName+"_fd_a0_r "+lineName+"_chi2 "+lineName+"_ndof"
		headerPV=" "+lineName+"_a0 "+lineName+"_flux "+lineName+"_fluxErr "+lineName+"_sigma "+lineName+"_sigmaErr "+lineName+"_continu "+lineName+"_continuErr "+lineName+"_EW "+lineName+"_share "+lineName+"_shareErr "+lineName+"_fd_a0_l "+lineName+"_fd_a0_r "+lineName+"_chi2 "+lineName+"_ndof"
		outPutNF=n.array([a0, self.dV,self.dV,self.dV, self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV])
		outPutNF_PV=n.array([a0, self.dV,self.dV,self.dV, self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV])
		modNF=n.array([self.dV,self.dV])
		if continuumSide=="left":
			domainLine=(wl>a0-fitWidth)&(wl<a0+fitWidth)
			domainCont=(wl>a0-DLC-fitWidth)&(wl<a0-fitWidth)
			if a0<wl.max()-DLC and a0>wl.min()+fitWidth and a0<wl.max()-fitWidth and len(domainLine.nonzero()[0])>2 and len(domainCont.nonzero()[0])>2 :
				continu=n.median(spec1d[domainCont])
				continuErr=n.median(err1d[domainCont])
				if model=="gaussian":
					flMod=lambda aa,sigma,F0,a0,continu : self.gaussianLine(aa,sigma,F0,a0,continu)
					p0=n.array([p0_sigma,p0_flux,a0,continu])
				if model=="lorentz":
					flMod=lambda aa,sigma,F0,a0,continu : self.lorentzLine(aa,sigma,F0,a0,continu)
					p0=n.array([p0_sigma,p0_flux, a0,continu])
				if model=="pseudoVoigt":
					flMod=lambda aa,sigma,F0,sh,a0,continu : self.pseudoVoigtLine(aa,sigma,F0,a0,continu,sh)
					p0=n.array([p0_sigma,p0_flux,p0_share,a0,continu])
				interp=interp1d(wl,spec1d)
				fd_a0_r=interp(a0+0.2)
				fd_a0_l=interp(a0-0.2)
				if fd_a0_r>continu and fd_a0_l>continu :
					out = curve_fit(flMod, wl[domainLine], spec1d[domainLine], p0=p0,sigma=err1d[domainLine],maxfev=1000000000, gtol=1.49012e-8)
					if out[1].__class__==n.ndarray and ( model=="gaussian" or model=="lorentz") : 
						model1=flMod(wl[domainLine],out[0][0],out[0][1],out[0][2],out[0][3])
						var=err1d[domainLine]
						chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
						ndof=len(var)
						sigma=out[0][0]
						sigmaErr=out[1][0][0]**0.5
						flux=out[0][1]
						fluxErr=out[1][1][1]**0.5
						a0=out[0][2]
						a0_err=out[1][2][2]**0.5
						continu=out[0][3]
						continuErr=out[1][3][3]**0.5
						EW=flux/continu
						outPut=n.array([a0,flux,fluxErr,sigma,sigmaErr,continu,continuErr,EW,fd_a0_l,fd_a0_r,chi2,ndof ])
						mod=n.array([wl[domainLine],model1])
						return outPut,mod,header
					elif model=="gaussian" or model=="lorentz" :
						return n.array([a0,self.dV,self.dV,self.dV, self.dV,continu,continuErr,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV ]),modNF,header
					elif out[1].__class__==n.ndarray and model=="pseudoVoigt" : 
						model1=flMod(wl[domainLine],out[0][0],out[0][1],out[0][2],out[0][3],out[0][4])
						var=err1d[domainLine]
						chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
						ndof=len(var)
						sigma=out[0][0]
						sigmaErr=out[1][0][0]**0.5
						flux=out[0][1]
						fluxErr=out[1][1][1]**0.5
						share=out[0][2]
						shareErr=out[1][2][2]**0.5
						a0=out[0][3]
						a0_err=out[1][3][3]**0.5
						continu=out[0][4]
						continuErr=out[1][4][4]**0.5
						EW=flux/continu
						outPut=n.array([a0,flux,fluxErr,sigma,sigmaErr,continu,continuErr,EW,share,shareErr,fd_a0_l,fd_a0_r,chi2,ndof])
						mod=n.array([wl[domainLine],model1])
						return outPut,mod,headerPV
					else :
						return n.array([a0,self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,self.dV,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,headerPV
				else :
					if  model=="gaussian" or model=="lorentz" :
						return n.array([a0,self.dV,self.dV,self.dV, self.dV,continu,continuErr,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV ]),modNF,header
					if model=="pseudoVoigt" :
						return n.array([a0,self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,self.dV,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,headerPV
			else :
				print "not enough space to fit the line"
				if  model=="gaussian" or model=="lorentz" :
					return outPutNF,modNF,header
				if model=="pseudoVoigt" :
					return outPutNF_PV,modNF,headerPV

		elif continuumSide=="right" :
			domainLine=(wl>a0-fitWidth)&(wl<a0+fitWidth)
			domainCont=(wl>a0+fitWidth)&(wl<a0+DLC+fitWidth)
			if a0<wl.max()-DLC and a0>wl.min()+fitWidth and a0<wl.max()-fitWidth and len(domainLine.nonzero()[0])>2 and len(domainCont.nonzero()[0])>2 :
				continu=n.median(spec1d[domainCont])
				continuErr=n.median(err1d[domainCont])
				if model=="gaussian":
					flMod=lambda aa,sigma,F0,a0,continu : self.gaussianLine(aa,sigma,F0,a0,continu)
					p0=n.array([p0_sigma,p0_flux,a0,continu])
				if model=="lorentz":
					flMod=lambda aa,sigma,F0,a0,continu : self.lorentzLine(aa,sigma,F0,a0,continu)
					p0=n.array([p0_sigma,p0_flux, a0,continu])
				if model=="pseudoVoigt":
					flMod=lambda aa,sigma,F0,sh,a0,continu : self.pseudoVoigtLine(aa,sigma,F0,a0,continu,sh)
					p0=n.array([p0_sigma,p0_flux,p0_share,a0,continu])
				interp=interp1d(wl,spec1d)
				fd_a0_r=interp(a0+0.2)
				fd_a0_l=interp(a0-0.2)
				if fd_a0_r>continu and fd_a0_l>continu :
					out = curve_fit(flMod, wl[domainLine], spec1d[domainLine], p0=p0,sigma=err1d[domainLine],maxfev=1000000000, gtol=1.49012e-8)
					if out[1].__class__==n.ndarray and ( model=="gaussian" or model=="lorentz") : 
						model1=flMod(wl[domainLine],out[0][0],out[0][1],out[0][2],out[0][3])
						var=err1d[domainLine]
						chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
						ndof=len(var)
						sigma=out[0][0]
						sigmaErr=out[1][0][0]**0.5
						flux=out[0][1]
						fluxErr=out[1][1][1]**0.5
						a0=out[0][2]
						a0_err=out[1][2][2]**0.5
						continu=out[0][3]
						continuErr=out[1][3][3]**0.5
						EW=flux/continu
						outPut=n.array([a0,flux,fluxErr,sigma,sigmaErr,continu,continuErr,EW,fd_a0_l,fd_a0_r,chi2,ndof ])
						mod=n.array([wl[domainLine],model1])
						return outPut,mod,header
					elif model=="gaussian" or model=="lorentz" :
						return n.array([a0,self.dV,self.dV,self.dV, self.dV,continu,continuErr,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV ]),modNF,header
					elif out[1].__class__==n.ndarray and model=="pseudoVoigt" : 
						model1=flMod(wl[domainLine],out[0][0],out[0][1],out[0][2],out[0][3],out[0][4])
						var=err1d[domainLine]
						chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
						ndof=len(var)
						sigma=out[0][0]
						sigmaErr=out[1][0][0]**0.5
						flux=out[0][1]
						fluxErr=out[1][1][1]**0.5
						share=out[0][2]
						shareErr=out[1][2][2]**0.5
						a0=out[0][3]
						a0_err=out[1][3][3]**0.5
						continu=out[0][4]
						continuErr=out[1][4][4]**0.5
						EW=flux/continu
						outPut=n.array([a0,flux,fluxErr,sigma,sigmaErr,continu,continuErr,EW,share,shareErr,fd_a0_l,fd_a0_r,chi2,ndof])
						mod=n.array([wl[domainLine],model1])
						return outPut,mod,headerPV
					else :
						return n.array([a0,self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,self.dV,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,headerPV
				else :
					if  model=="gaussian" or model=="lorentz" :
						return n.array([a0,self.dV,self.dV,self.dV, self.dV,continu,continuErr,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV ]),modNF,header
					if model=="pseudoVoigt" :
						return n.array([a0,self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,self.dV,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,headerPV
			else :
				print "not enough space to fit the line"
				if  model=="gaussian" or model=="lorentz" :
					return outPutNF,modNF,header
				if model=="pseudoVoigt" :
					return outPutNF_PV,modNF,headerPV

	def fit_Line(self,wl,spec1d,err1d,a0,lineName="AL",fitWidth=20,DLC=20, p0_sigma=15.,p0_flux=8e-17,p0_share=0.5,continuumSide="left",model="gaussian"):
		"""
		fits a line profile to a spectrum around a fixed line position

		:param wl: wavelength (array, Angstrom)
		:param spec1d: flux observed in the broad band (array, f lambda)
		:param err1d: error on the flux observed in the broad band (array, f lambda)
		:param a0: expected position of the peak of the line in the observed frame (redshifted). a0 is not fitted, it is given.
		:param lineName: suffix characterizing the line in the headers of the output
		:param DLC: wavelength extent to fit the continuum around the line. (def: 230 Angstrom)
		:param p0_sigma: prior on the line width in A (def: 15 A)
		:param p0_flux: prior on the line flux in erg/cm2/s/A (def: 8e-17)
		:param p0_share: prior on the share of Gaussian and Lorentzian model. Only used if the line is fitted with a pseudoVoigt profile width (def: 0.5 no units)
		:param continuumSide: "left" = bluewards of the line or "right" = redwards of the line
		:param model: line model to be fitted : "gaussian", "lorentz" or "pseudoVoigt".

		Returns :
		 * array 1 with the parameters of the model
		 * array 2 with the model (wavelength, flux model)
		 * header corresponding to the array 1
		"""
		header=" "+lineName+"_a0 "+lineName+"_flux "+lineName+"_fluxErr "+lineName+"_sigma "+lineName+"_sigmaErr "+lineName+"_continu "+lineName+"_continuErr "+lineName+"_EW "+lineName+"_fd_a0_l "+lineName+"_fd_a0_r "+lineName+"_chi2 "+lineName+"_ndof"
		headerPV=" "+lineName+"_a0 "+lineName+"_flux "+lineName+"_fluxErr "+lineName+"_sigma "+lineName+"_sigmaErr "+lineName+"_continu "+lineName+"_continuErr "+lineName+"_EW "+lineName+"_share "+lineName+"_shareErr "+lineName+"_fd_a0_l "+lineName+"_fd_a0_r "+lineName+"_chi2 "+lineName+"_ndof"
		outPutNF=n.array([a0, self.dV,self.dV,self.dV, self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV])
		outPutNF_PV=n.array([a0, self.dV,self.dV,self.dV, self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV])
		modNF=n.array([self.dV,self.dV])
		if continuumSide=="left":
			domainLine=(wl>a0-fitWidth)&(wl<a0+fitWidth)
			domainCont=(wl>a0-DLC-fitWidth)&(wl<a0-fitWidth)
			if a0<wl.max()-DLC and a0>wl.min()+fitWidth and a0<wl.max()-fitWidth and len(domainLine.nonzero()[0])>2 and len(domainCont.nonzero()[0])>2 :
				continu=n.median(spec1d[domainCont])
				continuErr=n.median(err1d[domainCont])
				if model=="gaussian":
					flMod=lambda aa,sigma,F0 : self.gaussianLine(aa,sigma,F0,a0,continu)
					p0=n.array([p0_sigma,p0_flux])
				if model=="lorentz":
					flMod=lambda aa,sigma,F0 : self.lorentzLine(aa,sigma,F0,a0,continu)
					p0=n.array([p0_sigma,p0_flux])
				if model=="pseudoVoigt":
					flMod=lambda aa,sigma,F0,sh : self.pseudoVoigtLine(aa,sigma,F0,a0,continu,sh)
					p0=n.array([p0_sigma,p0_flux,p0_share])
				interp=interp1d(wl,spec1d)
				fd_a0_r=interp(a0+0.2)
				fd_a0_l=interp(a0-0.2)
				if fd_a0_r>continu and fd_a0_l>continu :
					out = curve_fit(flMod, wl[domainLine], spec1d[domainLine], p0=p0,sigma=err1d[domainLine],maxfev=1000000000, gtol=1.49012e-8)
					if out[1].__class__==n.ndarray and ( model=="gaussian" or model=="lorentz") : 
						model1=flMod(wl[domainLine],out[0][0],out[0][1])
						var=err1d[domainLine]
						chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
						ndof=len(var)
						sigma=out[0][0]
						sigmaErr=out[1][0][0]**0.5
						flux=out[0][1]
						fluxErr=out[1][1][1]**0.5
						EW=flux/continu
						outPut=n.array([a0,flux,fluxErr,sigma,sigmaErr,continu,continuErr,EW,fd_a0_l,fd_a0_r,chi2,ndof ])
						mod=n.array([wl[domainLine],model1])
						return outPut,mod,header
					elif model=="gaussian" or model=="lorentz" :
						return n.array([a0,self.dV,self.dV,self.dV, self.dV,continu,continuErr,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV ]),modNF,header
					elif out[1].__class__==n.ndarray and model=="pseudoVoigt" : 
						model1=flMod(wl[domainLine],out[0][0],out[0][1],out[0][2])
						var=err1d[domainLine]
						chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
						ndof=len(var)
						sigma=out[0][0]
						sigmaErr=out[1][0][0]**0.5
						flux=out[0][1]
						fluxErr=out[1][1][1]**0.5
						share=out[0][2]
						shareErr=out[1][2][2]**0.5
						EW=flux/continu
						outPut=n.array([a0,flux,fluxErr,sigma,sigmaErr,continu,continuErr,EW,share,shareErr,fd_a0_l,fd_a0_r,chi2,ndof])
						mod=n.array([wl[domainLine],model1])
						return outPut,mod,headerPV
					else :
						return n.array([a0,self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,self.dV,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,headerPV
				else :
					if  model=="gaussian" or model=="lorentz" :
						return n.array([a0,self.dV,self.dV,self.dV, self.dV,continu,continuErr,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV ]),modNF,header
					if model=="pseudoVoigt" :
						return n.array([a0,self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,self.dV,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,headerPV
			else :
				print "not enough space to fit the line"
				if  model=="gaussian" or model=="lorentz" :
					return outPutNF,modNF,header
				if model=="pseudoVoigt" :
					return outPutNF_PV,modNF,headerPV

		elif continuumSide=="right" :
			domainLine=(wl>a0-fitWidth)&(wl<a0+fitWidth)
			domainCont=(wl>a0+fitWidth)&(wl<a0+DLC+fitWidth)
			if a0<wl.max()-DLC and a0>wl.min()+fitWidth and a0<wl.max()-fitWidth and len(domainLine.nonzero()[0])>2 and len(domainCont.nonzero()[0])>2 :
				continu=n.median(spec1d[domainCont])
				continuErr=n.median(err1d[domainCont])
				if model=="gaussian":
					flMod=lambda aa,sigma,F0 : self.gaussianLine(aa,sigma,F0,a0,continu)
					p0=n.array([p0_sigma,p0_flux])
				if model=="lorentz":
					flMod=lambda aa,sigma,F0 : self.lorentzLine(aa,sigma,F0,a0,continu)
					p0=n.array([p0_sigma,p0_flux])
				if model=="pseudoVoigt":
					flMod=lambda aa,sigma,F0,sh : self.pseudoVoigtLine(aa,sigma,F0,a0,continu,sh)
					p0=n.array([p0_sigma,p0_flux,p0_share])
				interp=interp1d(wl,spec1d)
				fd_a0_r=interp(a0+0.2)
				fd_a0_l=interp(a0-0.2)
				if fd_a0_r>continu and fd_a0_l>continu :
					out = curve_fit(flMod, wl[domainLine], spec1d[domainLine], p0=p0,sigma=err1d[domainLine],maxfev=1000000000, gtol=1.49012e-8)
					if out[1].__class__==n.ndarray and ( model=="gaussian" or model=="lorentz") : 
						model1=flMod(wl[domainLine],out[0][0],out[0][1])
						var=err1d[domainLine]
						chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
						ndof=len(var)
						sigma=out[0][0]
						sigmaErr=out[1][0][0]**0.5
						flux=out[0][1]
						fluxErr=out[1][1][1]**0.5
						EW=flux/continu
						outPut=n.array([a0,flux,fluxErr,sigma,sigmaErr,continu,continuErr,EW,fd_a0_l,fd_a0_r,chi2,ndof ])
						mod=n.array([wl[domainLine],model1])
						return outPut,mod,header
					elif model=="gaussian" or model=="lorentz" :
						return n.array([a0,self.dV,self.dV,self.dV, self.dV,continu,continuErr,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV ]),modNF,header
					elif out[1].__class__==n.ndarray and model=="pseudoVoigt" : 
						model1=flMod(wl[domainLine],out[0][0],out[0][1],out[0][2])
						var=err1d[domainLine]
						chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
						ndof=len(var)
						sigma=out[0][0]
						sigmaErr=out[1][0][0]**0.5
						flux=out[0][1]
						fluxErr=out[1][1][1]**0.5
						share=out[0][2]
						shareErr=out[1][2][2]**0.5
						EW=flux/continu
						outPut=n.array([a0,flux,fluxErr,sigma,sigmaErr,continu,continuErr,EW,share,shareErr,fd_a0_l,fd_a0_r,chi2,ndof])
						mod=n.array([wl[domainLine],model1])
						return outPut,mod,headerPV
					else :
						return n.array([a0,self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,self.dV,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,headerPV
				else :
					if  model=="gaussian" or model=="lorentz" :
						return n.array([a0,self.dV,self.dV,self.dV, self.dV,continu,continuErr,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV ]),modNF,header
					if model=="pseudoVoigt" :
						return n.array([a0,self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,self.dV,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,headerPV
			else :
				print "not enough space to fit the line"
				if  model=="gaussian" or model=="lorentz" :
					return outPutNF,modNF,header
				if model=="pseudoVoigt" :
					return outPutNF_PV,modNF,headerPV

	def fit_Line_OIIdoublet(self,wl,spec1d,err1d,a0=3726.0321735398957,lineName="OII",fitWidth=20,DLC=20,p0_sigma=4.,p0_flux=1e-16,p0_share=0.58,model="gaussian"):
		"""
		fits the [OII] doublet line profile

		:param wl: wavelength (array, Angstrom)
		:param spec1d: flux observed in the broad band (array, f lambda)
		:param err1d: error on the flux observed in the broad band (array, f lambda)
		:param a0: expected position of the peak of the line in the observed frame (redshifted). 2 positions given.
		:param lineName: suffix characterizing the line in the headers of the output
		:param DLC: wavelength extent to fit the continuum around the line. (def: 230 Angstrom)
		:param p0_sigma: prior on the line width in A (def: 15 A)
		:param p0_flux: prior on the line flux in erg/cm2/s/A (def: 8e-17)
		:param p0_share: prior on the share between the two [OII] lines. (def: 0.58)
		:param continuumSide: "left" = bluewards of the line or "right" = redwards of the line
		:param model: line model to be fitted : "gaussian", "lorentz" or "pseudoVoigt"

		Returns :
		 * array 1 with the parameters of the model
		 * array 2 with the model (wavelength, flux model)
		 * header corresponding to the array 1
		"""		
		header=" "+lineName+"_a0a "+lineName+"_a0b "+lineName+"_flux "+lineName+"_fluxErr "+lineName+"_sigma "+lineName+"_sigmaErr "+lineName+"_continu "+lineName+"_continuErr "+lineName+"_EW "+lineName+"_share "+lineName+"_shareErr "+lineName+"_fd_a0_l "+lineName+"_fd_a0_r "+lineName+"_chi2 "+lineName+"_ndof" 
		outPutNF=n.array([a0[0], a0[1], self.dV,self.dV, self.dV,self.dV, self.dV, self.dV,self.dV, self.dV,self.dV, self.dV,self.dV,self.dV,self.dV])
		modNF=n.array([self.dV,self.dV])
		domainLine=(wl>a0[0]-fitWidth)&(wl<a0[1]+fitWidth)
		domainCont=(wl>a0[0]-DLC-fitWidth)&(wl<a0[0]-fitWidth)
		if a0[0]<wl.max()-DLC and len(domainLine.nonzero()[0])>2 and len(domainCont.nonzero()[0])>2 :
			continu=n.median(spec1d[domainCont])
			continuErr=n.median(err1d[domainCont])
			if model=="gaussian":
				flMod=lambda aa,sigma,F0,sh :continu+ self.gaussianLineNC(aa,sigma,(1-sh)*F0,a0[0])+self.gaussianLineNC(aa,sigma,sh*F0,a0[1])
			if model=="lorentz":
				flMod=lambda aa,sigma,F0,sh : self.lorentzLine(aa,sigma,(1-sh)*F0,a0[0],continu/2.)+self.lorentzLine(aa,sigma,sh*F0,a0[1],continu/2.)
			index=n.searchsorted(wl,a0[1])
			fd_a0_r=spec1d[index]
			fd_a0_l=spec1d[index]
			index=n.searchsorted(wl,a0[0])
			if fd_a0_r>continu or fd_a0_l>continu :
				out = curve_fit(flMod, wl[domainLine], spec1d[domainLine], p0=n.array([p0_sigma,p0_flux,p0_share]),sigma=err1d[domainLine],maxfev=1000000000, gtol=1.49012e-8)
				if out[1].__class__==n.ndarray : 
					model1=flMod(wl[domainLine],out[0][0],out[0][1],out[0][2])
					var=err1d[domainLine]
					chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
					ndof=len(var)
					sigma=out[0][0]
					sigmaErr=out[1][0][0]**0.5
					flux=out[0][1]
					fluxErr=out[1][1][1]**0.5
					share=out[0][2]
					shareErr=out[1][2][2]**0.5
					EW=flux/continu
					outPut=n.array([a0[0],a0[1],flux,fluxErr,sigma,sigmaErr,continu,continuErr,EW,share,shareErr,fd_a0_l,fd_a0_r,chi2,ndof])
					mod=n.array([wl[domainLine],model1])
					return outPut,mod,header
				else :
					return n.array([a0[0],a0[1],self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,self.dV,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,header
			else :
				return n.array([a0[0],a0[1],self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,self.dV,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,header
		else :
			print "not enough space to fit the line"
			return outPutNF,modNF,header

	def fit_Line_OIIdoublet_position(self,wl,spec1d,err1d,a0=3726.0321,lineName="O2_3728",fitWidth=20,DLC=20,p0_sigma=4.,p0_flux=1e-16,p0_share=0.58,model="gaussian"):
		"""
		fits the [OII] doublet line profile

		:param wl: wavelength (array, Angstrom)
		:param spec1d: flux observed in the broad band (array, f lambda)
		:param err1d: error on the flux observed in the broad band (array, f lambda)
		:param a0: expected position of the peak of the line in the observed frame (redshifted). 2 positions given.
		:param lineName: suffix characterizing the line in the headers of the output
		:param DLC: wavelength extent to fit the continuum around the line. (def: 230 Angstrom)
		:param p0_sigma: prior on the line width in A (def: 15 A)
		:param p0_flux: prior on the line flux in erg/cm2/s/A (def: 8e-17)
		:param p0_share: prior on the share between the two [OII] lines. (def: 0.58)
		:param continuumSide: "left" = bluewards of the line or "right" = redwards of the line
		:param model: line model to be fitted : "gaussian", "lorentz" or "pseudoVoigt"

		Returns :
		 * array 1 with the parameters of the model
		 * array 2 with the model (wavelength, flux model)
		 * header corresponding to the array 1
		"""		
		header=" "+lineName+"_a0a "+lineName+"_a0b "+lineName+"_flux "+lineName+"_fluxErr "+lineName+"_sigma "+lineName+"_sigmaErr "+lineName+"_continu "+lineName+"_continuErr "+lineName+"_EW "+lineName+"_share "+lineName+"_shareErr "+lineName+"_fd_a0_l "+lineName+"_fd_a0_r "+lineName+"_chi2 "+lineName+"_ndof" 
		outPutNF=n.array([a0, a0+2.782374, self.dV,self.dV, self.dV,self.dV, self.dV, self.dV,self.dV, self.dV,self.dV, self.dV,self.dV,self.dV,self.dV])
		modNF=n.array([self.dV,self.dV])
		domainLine=(wl>a0-fitWidth)&(wl<a0+2.782374+fitWidth/2.)
		domainCont=(wl>a0-fitWidth-DLC)&(wl<a0-fitWidth)
		if a0<wl.max()-DLC and len(domainLine.nonzero()[0])>2 and len(domainCont.nonzero()[0])>2 :
			continu=n.median(spec1d[domainCont])
			continuErr=n.median(err1d[domainCont])
			if model=="gaussian":
				flMod=lambda aa,sigma,F0,sh,a0,continu :continu+ self.gaussianLineNC(aa,sigma,(1-sh)*F0,a0)+self.gaussianLineNC(aa,sigma,sh*F0,a0+2.782374)
				p0=n.array([p0_sigma,p0_flux,p0_share,a0,continu])
			index=n.searchsorted(wl,a0+2.782374)
			fd_a0_r=spec1d[index]
			index=n.searchsorted(wl,a0)
			fd_a0_l=spec1d[index]
			if fd_a0_r>continu or fd_a0_l>continu :
				out = curve_fit(flMod, wl[domainLine], spec1d[domainLine], p0=n.array([p0_sigma,p0_flux,p0_share,a0,continu]),sigma=err1d[domainLine],maxfev=1000000000, gtol=1.49012e-8)
				if out[1].__class__==n.ndarray : 
					model1=flMod(wl[domainLine],out[0][0],out[0][1],out[0][2],out[0][3],out[0][4])
					var=err1d[domainLine]
					chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
					ndof=len(var)
					sigma=out[0][0]
					sigmaErr=out[1][0][0]**0.5
					flux=out[0][1]
					fluxErr=out[1][1][1]**0.5
					share=out[0][2]
					shareErr=out[1][2][2]**0.5
					a0=out[0][3]
					a0_err=out[1][3][3]**0.5
					continu=out[0][4]
					continuErr=out[1][4][4]**0.5
					EW=flux/continu
					outPut=n.array([a0,a0+2.782374,flux,fluxErr,sigma,sigmaErr,continu,continuErr,EW,share,shareErr,fd_a0_l,fd_a0_r,chi2,ndof])
					mod=n.array([wl[domainLine],model1])
					return outPut,mod,header
				else :
					return n.array([a0,a0+2.782374,self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,self.dV,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,header
			else :
				return n.array([a0,a0+2.782374,self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,self.dV,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,header
		else :
			print "not enough space to fit the line"
			return outPutNF,modNF,header

	def fit_Line_OIIdoublet_position_C0noise(self,wl,spec1d,err1d,a0=3726.0321,lineName="O2_3728",fitWidth=20,DLC=20,p0_sigma=4.,p0_flux=1e-16,p0_share=0.58,model="gaussian"):
		"""
		fits the [OII] doublet line profile

		:param wl: wavelength (array, Angstrom)
		:param spec1d: flux observed in the broad band (array, f lambda)
		:param err1d: error on the flux observed in the broad band (array, f lambda)
		:param a0: expected position of the peak of the line in the observed frame (redshifted). 2 positions given.
		:param lineName: suffix characterizing the line in the headers of the output
		:param DLC: wavelength extent to fit the continuum around the line. (def: 230 Angstrom)
		:param p0_sigma: prior on the line width in A (def: 15 A)
		:param p0_flux: prior on the line flux in erg/cm2/s/A (def: 8e-17)
		:param p0_share: prior on the share between the two [OII] lines. (def: 0.58)
		:param continuumSide: "left" = bluewards of the line or "right" = redwards of the line
		:param model: line model to be fitted : "gaussian", "lorentz" or "pseudoVoigt"

		Returns :
		 * array 1 with the parameters of the model
		 * array 2 with the model (wavelength, flux model)
		 * header corresponding to the array 1
		"""		
		header=" "+lineName+"_a0a "+lineName+"_a0b "+lineName+"_flux "+lineName+"_fluxErr "+lineName+"_sigma "+lineName+"_sigmaErr "+lineName+"_continu "+lineName+"_continuErr "+lineName+"_EW "+lineName+"_share "+lineName+"_shareErr "+lineName+"_fd_a0_l "+lineName+"_fd_a0_r "+lineName+"_chi2 "+lineName+"_ndof" 
		outPutNF=n.array([a0, a0+2.782374, self.dV,self.dV, self.dV,self.dV, self.dV, self.dV,self.dV, self.dV,self.dV, self.dV,self.dV,self.dV,self.dV])
		modNF=n.array([self.dV,self.dV])
		domainLine=(wl>a0-fitWidth)&(wl<a0+2.782374+fitWidth/2.)
		domainCont=(wl>a0-fitWidth-DLC)&(wl<a0-fitWidth)
		if a0<wl.max()-DLC and len(domainLine.nonzero()[0])>2 and len(domainCont.nonzero()[0])>2 :
			continu=n.median(spec1d[domainCont])
			continuErr=n.median(err1d[domainCont])
			if model=="gaussian":
				flMod=lambda aa,sigma,F0,sh,a0,continu :continu+ self.gaussianLineNC(aa,sigma,(1-sh)*F0,a0)+self.gaussianLineNC(aa,sigma,sh*F0,a0+2.782374)
				p0=n.array([p0_sigma,p0_flux,p0_share,a0,continu])
			index=n.searchsorted(wl,a0+2.782374)
			fd_a0_r=spec1d[index]
			index=n.searchsorted(wl,a0)
			fd_a0_l=spec1d[index]
			if fd_a0_r>continu or fd_a0_l>continu :
				out = curve_fit(flMod, wl[domainLine], spec1d[domainLine], p0=n.array([p0_sigma,p0_flux,p0_share,a0,continu]),sigma=continu*n.ones_like(err1d[domainLine]),maxfev=1000000000, gtol=1.49012e-8)
				if out[1].__class__==n.ndarray : 
					model1=flMod(wl[domainLine],out[0][0],out[0][1],out[0][2],out[0][3],out[0][4])
					var=continu*n.ones_like(err1d[domainLine])
					chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
					ndof=len(var)
					sigma=out[0][0]
					sigmaErr=out[1][0][0]**0.5
					flux=out[0][1]
					fluxErr=out[1][1][1]**0.5
					share=out[0][2]
					shareErr=out[1][2][2]**0.5
					a0=out[0][3]
					a0_err=out[1][3][3]**0.5
					continu=out[0][4]
					continuErr=out[1][4][4]**0.5
					EW=flux/continu
					outPut=n.array([a0,a0+2.782374,flux,fluxErr,sigma,sigmaErr,continu,continuErr,EW,share,shareErr,fd_a0_l,fd_a0_r,chi2,ndof])
					mod=n.array([wl[domainLine],model1])
					return outPut,mod,header
				else :
					return n.array([a0,a0+2.782374,self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,self.dV,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,header
			else :
				return n.array([a0,a0+2.782374,self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,self.dV,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,header
		else :
			print "not enough space to fit the line"
			return outPutNF,modNF,header

	def fit_recLine(self,wl,spec1d,err1d,a0,lineName="AL",fitWidth=20,DLC=20,p0_sigma=5.,p0_flux=5e-17,continuumSide="left"):
		"""
		fits a recombination line profile : emission and absorption modeled by Gaussians. Only for high SNR spectra.

		:param wl: wavelength (array, Angstrom)
		:param spec1d: flux observed in the broad band (array, f lambda)
		:param err1d: error on the flux observed in the broad band (array, f lambda)
		:param a0: expected position of the peak of the line in the observed frame (redshifted)
		:param lineName: suffix characterizing the line in the headers of the output
		:param DLC: wavelength extent to fit the continuum around the line. (def: 230 Angstrom)
		:param p0_sigma: prior on the line width in A (def: 15 A)
		:param p0_flux: prior on the line flux in erg/cm2/s/A (def: 8e-17)
		:param continuumSide: "left" = bluewards of the line or "right" = redwards of the line

		Returns :
		 * array 1 with the parameters of the model
		 * array 2 with the model (wavelength, flux model)
		 * header corresponding to the array 1
		"""		
		header=" "+lineName+"_a0 "+lineName+"_flux "+lineName+"_fluxErr "+lineName+"_sigma "+lineName+"_sigmaErr "+lineName+"_continu "+lineName+"_continuErr "+lineName+"_EW "+lineName+"_fd_a0_l "+lineName+"_fd_a0_r "+lineName+"_chi2 "+lineName+"_ndof"
		headerPV=" "+lineName+"_a0 "+lineName+"_flux "+lineName+"_fluxErr "+lineName+"_sigma "+lineName+"_sigmaErr "+lineName+"_continu "+lineName+"_continuErr "+lineName+"_EW "+lineName+"_share "+lineName+"_shareErr_"+" fd_a0_l "+lineName+"_fd_a0_r "+lineName+"_chi2 "+lineName+"_ndof"
		outPutNF=n.array([a0, self.dV,self.dV,self.dV, self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV,self.dV])
		modNF=n.array([self.dV,self.dV])
		if continuumSide=="left":
			domainLine=(wl>a0-fitWidth)&(wl<a0+fitWidth)
			domainCont=(wl>a0-DLC)&(wl<a0-fitWidth)
			if a0<wl.max()-DLC and a0>wl.min()+fitWidth and a0<wl.max()-fitWidth and len(domainLine.nonzero()[0])>2 and len(domainCont.nonzero()[0])>2 :
				continu=n.median(spec1d[domainCont])
				continuErr=n.median(err1d[domainCont])
				# model with absorption
				flMod=lambda aa,sigma,F0,sigmaL,F0L,a0L,sigmaR,F0R,a0R : continu + self.gaussianLineNC(aa,sigma,F0,a0) - self.gaussianLineNC(aa,sigmaL,F0L,a0L) - self.gaussianLineNC(aa,sigmaR,F0R,a0R)
				p0=n.array([p0_sigma,p0_flux,p0_sigma/2.,p0_flux/5.,a0-5, p0_sigma/2.,p0_flux/5.,a0-5])
				interp=interp1d(wl,spec1d)
				fd_a0_r=interp(a0+0.2)
				fd_a0_l=interp(a0-0.2)
				if fd_a0_r>continu and fd_a0_l>continu :
					out = curve_fit(flMod, wl[domainLine], spec1d[domainLine], p0=p0,sigma=err1d[domainLine],maxfev=1000000000, gtol=1.49012e-8)
					if out[1].__class__==n.ndarray : 
						model1=flMod(wl[domainLine],out[0][0],out[0][1],out[0][2],out[0][3],out[0][4],out[0][5],out[0][6],out[0][7])
						var=err1d[domainLine]
						chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
						ndof=len(var)-8
						sigma=out[0][0]
						sigmaErr=out[1][0][0]**0.5
						flux=out[0][1]
						fluxErr=out[1][1][1]**0.5
						EW=flux/continu
						outPut=n.array([a0,flux,fluxErr,sigma,sigmaErr,continu,continuErr,EW,fd_a0_l,fd_a0_r,chi2,ndof ])
						mod=n.array([wl[domainLine],model1])
						return outPut,mod,header
					else :
						return n.array([a0,self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,headerPV
				else :
					return n.array([a0,self.dV,self.dV,self.dV, self.dV,continu,continuErr,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV ]),modNF,header
			else :
				print "not enough space to fit the line"
				return outPutNF,modNF,header

		elif continuumSide=="right" :
			domainLine=(wl>a0-fitWidth)&(wl<a0+fitWidth)
			domainCont=(wl>a0+fitWidth)&(wl<a0+DLC)
			if a0<wl.max()-DLC and a0>wl.min()+fitWidth and a0<wl.max()-fitWidth and len(domainLine.nonzero()[0])>2 and len(domainCont.nonzero()[0])>2 :
				continu=n.median(spec1d[domainCont])
				continuErr=n.median(err1d[domainCont])
				# model with absorption
				flMod=lambda aa,sigma,F0,sigmaL,F0L,a0L,sigmaR,F0R,a0R : continu + self.gaussianLineNC(aa,sigma,F0,a0) - self.gaussianLineNC(aa,sigmaL,F0L,a0L) - self.gaussianLineNC(aa,sigmaR,F0R,a0R)
				p0=n.array([p0_sigma,p0_flux,p0_sigma/2.,p0_flux/5.,a0-5, p0_sigma/2.,p0_flux/5.,a0-5])
				interp=interp1d(wl,spec1d)
				fd_a0_r=interp(a0+0.2)
				fd_a0_l=interp(a0-0.2)
				if fd_a0_r>continu and fd_a0_l>continu :
					out = curve_fit(flMod, wl[domainLine], spec1d[domainLine], p0=p0,sigma=err1d[domainLine],maxfev=1000000000, gtol=1.49012e-8)
					if out[1].__class__==n.ndarray : 
						model1=flMod(wl[domainLine],out[0][0],out[0][1],out[0][2],out[0][3],out[0][4],out[0][5],out[0][6],out[0][7])
						var=err1d[domainLine]
						chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
						ndof=len(var)-8
						sigma=out[0][0]
						sigmaErr=out[1][0][0]**0.5
						flux=out[0][1]
						fluxErr=out[1][1][1]**0.5
						EW=flux/continu
						outPut=n.array([a0,flux,fluxErr,sigma,sigmaErr,continu,continuErr,EW,fd_a0_l,fd_a0_r,chi2,ndof ])
						mod=n.array([wl[domainLine],model1])
						return outPut,mod,header
					else :
						return n.array([a0,self.dV,self.dV,self.dV,self.dV,continu,continuErr,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV]),modNF,headerPV
				else :
					return n.array([a0,self.dV,self.dV,self.dV, self.dV,continu,continuErr,self.dV,fd_a0_l,fd_a0_r,self.dV,self.dV ]),modNF,header
			else :
				print "not enough space to fit the line"
				return outPutNF,modNF,header


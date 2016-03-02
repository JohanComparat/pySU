"""
.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

General purpose:
................

The class ModelSpectraStacks is dedicated to modelling and extracting information from stacks of spectra.

*Imports*::

	import matplotlib
	matplotlib.use('pdf')
	import matplotlib.pyplot as p
	import os 
	import astropy.cosmology as co
	cosmo=co.FlatLambdaCDM(H0=70,Om0=0.3)
	import astropy.units as u
	import astropy.io.fits as fits
	import numpy as n
	from scipy.optimize import curve_fit
	from scipy.interpolate import interp1d
	from scipy.stats import scoreatpercentile
	import astropy.io.fits as fits
	from lineListAir import *
	import LineFittingLibrary as lineFit


"""
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p
import os 
import astropy.cosmology as co
cosmo=co.Planck15 #co.FlatLambdaCDM(H0=70,Om0=0.3)
import astropy.units as u
import astropy.io.fits as fits
import numpy as n
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.stats import scoreatpercentile

import astropy.io.fits as fits

from lineListAir import *
allLinesList = n.array([ [Ne3,Ne3_3869,"Ne3_3869","left"], [O3,O3_4363,"O3_4363","right"], [O3,O3_4960,"O3_4960","left"], [O3,O3_5007,"O3_5007","right"], [N2,N2_6549,"N2_6549","left"], [N2,N2_6585,"N2_6585","right"], [S2,S2_6718,"S2_6718","left"], [S2,S2_6732,"S2_6732","right"], [Ar3,Ar3_7137,"Ar3_7137","left"], [H1,H1_1216,"H1_1216","right"], [H1,H1_3970,"H1_3970","right"], [H1,H1_4102,"H1_4102","right"], [H1,H1_4341,"H1_4341","right"], [H1,H1_4862,"H1_4862","left"], [H1,H1_6564,"H1_6564","left"]])

doubletList = n.array([[O2_3727,"O2_3727",O2_3729,"O2_3729",O2_mean]])

# import the fitting routines
import LineFittingLibrary as lineFit

O2a=3727.092 
O2b=3729.875 
O2=(O2a+O2b)/2.
Hg=4102.892
Hd=4341.684
Hb=4862.683
O3a=4960.295
O3b=5008.240
Ha=6564.61

fnu = lambda mAB : 10**(-(mAB+48.6)/2.5) # erg/cm2/s/Hz
flambda= lambda mAB, ll : 10**10 * c*1000 * fnu(mAB) / ll**2. # erg/cm2/s/A

kla=lambda ll :2.659 *(-2.156+1.509/ll-0.198/ll**2+0.011/ll**3 ) + 4.05
klb=lambda ll :2.659 *(-1.857+1.040/ll)+4.05

def kl(ll):
	"""Calzetti extinction law"""
	if ll>6300:
		return klb(ll)
	if ll<=6300:
		return kla(ll)

klO2=kl(O2)
klO3=kl(O3b)
klHb=kl(Hb)

H1=pn.RecAtom('H',1) # Hydrogen Balmer series

bdc0_ref=H1.getEmissivity(1e4, 1e2, lev_i = 4, lev_j = 2) / H1.getEmissivity(1e4, 1e2, lev_i = 3, lev_j = 2)
bdc1_ref=H1.getEmissivity(1e4, 1e2, lev_i = 4, lev_j = 2) / H1.getEmissivity(1e4, 1e2, lev_i = 5, lev_j = 2)
bdc2_ref=H1.getEmissivity(1e4, 1e2, lev_i = 4, lev_j = 2) / H1.getEmissivity(1e4, 1e2, lev_i = 6, lev_j = 2)
bdc3_ref=H1.getEmissivity(1e4, 1e2, lev_i = 4, lev_j = 2) / H1.getEmissivity(1e4, 1e2, lev_i = 7, lev_j = 2)
bdc4_ref=H1.getEmissivity(1e4, 1e2, lev_i = 4, lev_j = 2) / H1.getEmissivity(1e4, 1e2, lev_i = 8, lev_j = 2)
bdc5_ref=H1.getEmissivity(1e4, 1e2, lev_i = 4, lev_j = 2) / H1.getEmissivity(1e4, 1e2, lev_i = 9, lev_j = 2)

bdc23_ref=H1.getEmissivity(1e4, 1e2, lev_i = 5, lev_j = 2)/H1.getEmissivity(1e4, 1e2, lev_i = 6, lev_j = 2)


class ModelSpectraStacks:
	"""
	This class fits the emission lines on the continuum-subtracted stack.

	:param stack_file: fits file generated with a LF in a luminosity bin.
	:param cosmo: cosmology class from astropy
	:param firefly_min_wavelength: minimum wavelength considered by firefly (default : 1000) 
	:param firefly_max_wavelength: minimum wavelength considered by firefly (default : 7500)
	:param dV: default value that hold the place (default : -9999.99) 
	:param N_spectra_limitFraction: If the stack was made with N spectra. N_spectra_limitFraction selects the points that have were computed using more thant N_spectra_limitFraction * N spectra. (default : 0.8)
	"""
	def __init__(self, stack_file, mode="MILES", cosmo=cosmo, firefly_min_wavelength= 1000., firefly_max_wavelength=7500., dV=-9999.99, N_spectra_limitFraction=0.8):
		self.stack_file = stack_file
		self.mode = mode
		if self.mode=="MILES":
			self.stack_model_file = self.stack_file[:-5]+"-SPM-MILES.fits"
		if self.mode=="STELIB":
			self.stack_model_file = self.stack_file[:-5]+"-SPM-STELIB.fits"

		self.cosmo = cosmo
		self.firefly_max_wavelength	= firefly_max_wavelength
		self.firefly_min_wavelength	= firefly_min_wavelength
		self.dV = dV
		self.side = ''
		self.N_spectra_limitFraction = N_spectra_limitFraction
		# define self.sphereCM, find redshift ...
		self.redshift = float(self.stack_file.split('-')[2].split('_')[0][1:])
		sphere=4*n.pi*( self.cosmo.luminosity_distance(self.redshift) )**2.
		self.sphereCM=sphere.to(u.cm**2)
		hdus = fits.open(self.stack_file)
		self.hdR = hdus[0].header
		self.hdu1 = hdus[1] # .data
		# loads the data :
		wlA,flA,flErrA = self.hdu1.data['wavelength'], self.hdu1.data['meanWeightedStack'], self.hdu1.data['jackknifStackErrors']
		self.selection = (flA>0) & (self.hdu1.data['NspectraPerPixel']  > float( self.stack_file.split('_')[-5]) * self.N_spectra_limitFraction )
		self.wl,self.fl,self.flErr = wlA[self.selection], flA[self.selection], flErrA[self.selection] 
		self.stack=interp1d(self.wl,self.fl)
		self.stackErr=interp1d(self.wl,self.flErr)
		# loads model :
		hdus = fits.open(self.stack_model_file)
		self.hdu2 = hdus[1] # .data
		self.wlModel,self.flModel = self.hdu2.data['wavelength'], self.hdu2.data['firefly_model']*10**(-17)
		self.model=interp1d(n.hstack((self.wlModel,[n.max(self.wlModel)+10,11000])), n.hstack(( self.flModel, [n.median(self.flModel[:-20]),n.median(self.flModel[:-20])] )) )
		# wavelength range common to the stack and the model :
		self.wlLineSpectrum  = n.arange(n.max([self.stack.x.min(),self.model.x.min()]), n.min([self.stack.x.max(),self.model.x.max()]), 0.5)[2:-1]
		self.flLineSpectrum=n.array([self.stack(xx)-self.model(xx) for xx in self.wlLineSpectrum])
		self.fl_frac_LineSpectrum=n.array([self.stack(xx)/self.model(xx) for xx in self.wlLineSpectrum])
		self.flErrLineSpectrum=self.stackErr(self.wlLineSpectrum)


	def interpolate_stack(self):
		"""
		Divides the measured stack in overlapping and non-overlapping parts with the model.
		"""
		self.stack=interp1d(self.wl,self.fl)
		self.stackErr=interp1d(self.wl,self.flErr)
		# bluer than model
		self.stBlue = (self.wl<=self.firefly_min_wavelength)
		# optical
		self.stOpt = (self.wl<self.firefly_max_wavelength)& (self.wl> self.firefly_min_wavelength) 
		# redder than model
		self.stRed = (self.wl>=self.firefly_max_wavelength) 
		if len(self.wl)<50 :
			print "no data, skips spectrum"
			return 0.
		if len(self.wl[self.stBlue])>0:
			self.contBlue=n.median(self.fl[self.stBlue])
			self.side='blue'
		if len(self.wl[self.stRed])>0:
			self.contRed=n.median(self.fl[self.stRed])
			self.side='red'
		if len(self.wl[self.stRed])>0 and len(self.wl[self.stBlue])>0:
			self.contRed=n.median(self.fl[self.stRed])
			self.contBlue=n.median(self.fl[self.stBlue])
			self.side='both'
		if len(self.wl[self.stRed])==0 and len(self.wl[self.stBlue])==0:
			self.side='none'

	def interpolate_model(self):
		"""
		Interpolates the model to an array with the same coverage as the stack.
		"""
		# overlap region with stack
		self.mdOK =(self.wlModel>n.min(self.wl))&(self.wlModel<n.max(self.wl)) 

		mdBlue=(self.wlModel<=n.min(self.wl)) # bluer part than data
		mdRed=(self.wlModel>=n.max(self.wl)) # redder part than data
		okRed=(self.wlModel>4650)&(self.wlModel<self.firefly_max_wavelength)

		# Correction model => stack
		CORRection=n.sum((self.wl[self.stOpt][1:]-self.wl[self.stOpt][:-1])* self.fl[self.stOpt][1:]) / n.sum((self.wlModel[ self.mdOK ][1:]-self.wlModel[ self.mdOK ][:-1])*   self.flModel [ self.mdOK ][1:])
		print "Correction", CORRection

		if self.side=='red':
			self.model=interp1d(n.hstack((self.wlModel[ self.mdOK ],n.arange(self.wlModel[ self.mdOK ].max()+0.5, stack.x.max(), 0.5))), n.hstack((  self.flModel [ self.mdOK ]*CORRection, n.ones_like(n.arange( self.wlModel[ self.mdOK ].max() + 0.5, stack.x.max(), 0.5))*contRed )) )
		elif self.side=='blue':
			self.model=interp1d(n.hstack((n.arange(stack.x.min(),self.wlModel[ self.mdOK ].min()-1., 0.5),self.wlModel[ self.mdOK ])),n.hstack(( n.ones_like(n.arange(stack.x.min() ,self.wlModel[ self.mdOK ].min() -1.,0.5))* contBlue,  self.flModel [ self.mdOK ]*CORRection )) )
		elif self.side=='both':
			x1=n.hstack((n.arange(stack.x.min(),self.wlModel[ self.mdOK ].min()-1., 0.5), self.wlModel[ self.mdOK ]))
			y1=n.hstack(( n.ones_like(n.arange(stack.x.min(),self.wlModel[ self.mdOK ].min()- 1.,0.5))*contBlue,  self.flModel [ self.mdOK ]*CORRection ))
			x2=n.hstack((x1,n.arange(self.wlModel[ self.mdOK ].max()+0.5,stack.x.max(),0.5)))
			y2=n.hstack((y1,n.ones_like(n.arange(self.wlModel[ self.mdOK ].max()+0.5, stack.x.max(), 0.5))*contRed )) 
			self.model=interp1d(x2,y2)

		elif self.side=='none':
			self.model=interp1d(self.wlModel[ self.mdOK ], self.flModel [ self.mdOK ])


	def subtract_continuum_model(self):
		"""
		Creates the continuum substracted spectrum: the 'line' spectrum.
		"""
		self.interpolate_stack()
		self.interpolate_model()
		# wavelength range common to the stack and the model :
		self.wlLineSpectrum  = n.arange(n.max([self.stack.x.min(),self.model.x.min()]), n.min([self.stack.x.max(),self.model.x.max()]), 0.5)[2:-1]
		print "range probed", self.wlLineSpectrum[0], self.wlLineSpectrum[-1], len( self.wlLineSpectrum)
		self.flLineSpectrum=n.array([self.stack(xx)-self.model(xx) for xx in self.wlLineSpectrum])
		self.flErrLineSpectrum=self.stackErr(self.wlLineSpectrum)

	def fit_lines_to_lineSpectrum(self):
		"""
		Fits the emission lines on the line spectrum.
		"""
		# interpolates the mean spectra.
		if self.stack_file.find('VVDS')>0 or self.stack_file.find('VIPERS')>0 :
			lfit  =  lineFit.LineFittingLibrary(fitWidth = 70.)
		if self.stack_file.find('DEEP2')>0 :
			lfit  =  lineFit.LineFittingLibrary(fitWidth = 40.)

		#self.subtract_continuum_model()
		data,h=[],[]

		dat_mean,mI,hI=lfit.fit_Line_OIIdoublet(self.wlLineSpectrum, self.flLineSpectrum, self.flErrLineSpectrum, a0= n.array([O2_3727,O2_3729]) , lineName="O2_3728", p0_sigma=1,model="gaussian")

		d_out=[]
		for kk in range(10):
			fluxRR = interp1d(self.wl, self.hdu1.data['jackknifeSpectra'].T[kk][self.selection])
			flLineSpectrumRR=n.array([fluxRR(xx)-self.model(xx) for xx in self.wlLineSpectrum])
			d1,mI,hI=lfit.fit_Line_OIIdoublet(self.wlLineSpectrum, flLineSpectrumRR, self.flErrLineSpectrum, a0= n.array([O2_3727,O2_3729]) , lineName="O2_3728", p0_sigma=1,model="gaussian")
			d_out.append(d1)

		d_out = n.array(d_out)
		err_out = n.std(d_out,axis=0)
		# assign error values :
		dat_mean[2] = err_out[2-1]
		dat_mean[4] = err_out[4-1]
		dat_mean[6] = err_out[6-1]
		data.append(dat_mean)
		h.append(hI)

		for li in allLinesList :
			# measure line properties from the mean weighted stack
			dat_mean,mI,hI=lfit.fit_Line(self.wl,self.fl,self.flErr, li[1], lineName=li[2], continuumSide=li[3], model="gaussian",p0_sigma=1)
			# measure its dispersion using the stacks
			d_out=[]
			for kk in range(len(self.hdu1.data['jackknifeSpectra'].T)):
				fluxRR = interp1d(self.wl, self.hdu1.data['jackknifeSpectra'].T[kk][self.selection])
				flLineSpectrumRR=n.array([fluxRR(xx)-self.model(xx) for xx in self.wlLineSpectrum])
				d1,mI,hI=lfit.fit_Line(self.wlLineSpectrum, flLineSpectrumRR, self.flErrLineSpectrum, li[1], lineName=li[2], continuumSide=li[3], model="gaussian",p0_sigma=1)
				d_out.append(d1)

			d_out = n.array(d_out)
			err_out = n.std(d_out,axis=0)
			# assign error values :
			dat_mean[2] = err_out[2-1]
			dat_mean[4] = err_out[4-1]
			dat_mean[6] = err_out[6-1]
			data.append(dat_mean)
			h.append(hI)

		heading="".join(h)
		out=n.hstack((data))
		out[n.isnan(out)]=n.ones_like(out[n.isnan(out)])*self.dV
		#output = n.array([ out ])
		#print "----------------", output.T[0], output.T[1], output
		colNames = heading.split()
		#col0 = fits.Column(name=colNames[0],format='D', array= output.T[0])
		#col1 = fits.Column(name=colNames[1],format='D', array= output.T[1])
		#cols = fits.ColDefs([col0, col1])
		#print colNames
		for ll in range(len(colNames)):
			self.hdR[colNames[ll]+"_nc"] = out.T[ll]
			#cols += fits.Column(name=colNames[ll],format='D', array= output.T[ll] )


	def compute_derived_quantities(self):
		"""
		Computes the different line ratios and converts to extinction :
		 * Balmer decrement and E(B-V) Correction using 4862 / 4341
		 * Balmer decrement and E(B-V) Correction using 4862 / 4102
		 * Balmer decrement and E(B-V) Correction using 4341 / 4102
		 * O32 : O3/O2
		 * intrinsic fluxes O3 4960 and 5007
		 * intrinsic fluxes O2
		 * intrinsic fluxes Hb
		 * SFR from [OII] 
		 * SFR from Hbeta
		 * compute R23
		 * compute R23 intrinsic
		 * 12 log(O/H) with Tremonti 04 estimator
		 * 12 log OH with O2(3728)/Hbeta
		 * 12 log OH with O3(4960+5007)/Hbeta
		 * 12 log OH with O3(5007)/Hbeta
		"""
		self.BDarray=n.array([0,0,0])
		if self.hdR['H1_4341_flux_nc']>0 and self.hdR['H1_4862_flux_nc']>0 and self.hdR['H1_4341_fluxErr_nc']>0 and self.hdR['H1_4862_fluxErr_nc']>0 :
			self.BDarray[0]=1
			# Balmer decrement : 4862 / 4341
			self.hdR['BD_4862_4341']=self.hdR['H1_4341_flux_nc']/ self.hdR['H1_4862_flux_nc']
			bdc1ErrFrac = ( (self.hdR['H1_4862_fluxErr_nc']/ self.hdR['H1_4862_flux_nc'])**2 + (self.hdR['H1_4341_fluxErr_nc']/ self.hdR['H1_4341_flux_nc'])**2. ) **0.5
			self.hdR['BD_4862_4341_err']= self.hdR['BD_4862_4341'] * bdc1ErrFrac
			# E(B-V) Correction using 4862 / 4341
			self.hdR['EBV_4862_4341'] = -5*n.log10(self.hdR['BD_4862_4341'] * bdc1_ref) / (2* (5.12 - 4.6))
			self.hdR['EBV_4862_4341_err']= -5 * bdc1ErrFrac * bdc1_ref/(2*(5.12-4.6)*n.log(10))
			# applied to emission lines using Calzetti's law
			self.hdR['EBV_4862_4341_CORRO2']=10**(0.4 * self.hdR['EBV_4862_4341'] *klO2)
			self.hdR['EBV_4862_4341_CORRO2_err']= self.hdR['EBV_4862_4341_err'] * n.log(10) * 0.4 * klO2 * self.hdR['EBV_4862_4341_CORRO2']
			self.hdR['EBV_4862_4341_CORRO3']= 10**(0.4 *self.hdR['EBV_4862_4341'] *klO3)
			self.hdR['EBV_4862_4341_CORRO3_err']= self.hdR['EBV_4862_4341_err'] * n.log(10) * 0.4 * klO3 * self.hdR['EBV_4862_4341_CORRO3']
			self.hdR['EBV_4862_4341_CORRHb']=10**(0.4 *klHb *self.hdR['EBV_4862_4341'] )
			self.hdR['EBV_4862_4341_CORRHb_err']= self.hdR['EBV_4862_4341_err'] * n.log(10) * 0.4 * klHb * self.hdR['EBV_4862_4341_CORRHb']
		else :
			self.hdR['BD_4862_4341'] = self.dV
			self.hdR['BD_4862_4341_err'] = self.dV
			self.hdR['EBV_4862_4341'] = self.dV
			self.hdR['EBV_4862_4341_err'] = self.dV
			self.hdR['EBV_4862_4341_CORRO2'] = self.dV
			self.hdR['EBV_4862_4341_CORRO2_err'] = self.dV
			self.hdR['EBV_4862_4341_CORRO3'] = self.dV
			self.hdR['EBV_4862_4341_CORRO3_err'] = self.dV
			self.hdR['EBV_4862_4341_CORRHb'] = self.dV
			self.hdR['EBV_4862_4341_CORRHb_err'] = self.dV

		if self.hdR['H1_4102_flux_nc']>0 and self.hdR['H1_4862_flux_nc']>0 and self.hdR['H1_4102_fluxErr_nc']>0 and self.hdR['H1_4862_fluxErr_nc']>0 :
			self.BDarray[1]=1
			# Balmer decrement : 4862 / 4102
			self.hdR['BD_4862_4102'] = self.hdR['H1_4102_flux_nc']/ self.hdR['H1_4862_flux_nc']
			bdc2ErrFrac = ( (self.hdR['H1_4862_fluxErr_nc']/ self.hdR['H1_4862_flux_nc'] )**2 + (self.hdR['H1_4102_fluxErr_nc']/self.hdR['H1_4102_flux_nc'])**2. ) **0.5
			self.hdR['BD_4862_4102_err'] = self.hdR['BD_4862_4102']* bdc2ErrFrac
			# E(B-V) Correction using 4862 / 4341
			self.hdR['EBV_4862_4102'] = -5*n.log10( self.hdR['BD_4862_4102'] * bdc2_ref )/( 2*(5.39-4.6))
			self.hdR['EBV_4862_4102_err'] = -5 * bdc2ErrFrac * bdc2_ref /(2*( 5.39 - 4.6)*n.log(10))
			# applied to emission lines using Calzetti's law
			self.hdR['EBV_4862_4102_CORRO2']=10**(0.4 *self.hdR['EBV_4862_4102'] *klO2)
			self.hdR['EBV_4862_4102_CORRO2_err']= self.hdR['EBV_4862_4102_err'] * n.log(10) * 0.4 * klO2 * self.hdR['EBV_4862_4102_CORRO2']
			self.hdR['EBV_4862_4102_CORRO3']= 10**(0.4 *self.hdR['EBV_4862_4102'] *klO3)
			self.hdR['EBV_4862_4102_CORRO3_err']= self.hdR['EBV_4862_4102_err'] * n.log(10) * 0.4 * klO3 * self.hdR['EBV_4862_4102_CORRO3']
			self.hdR['EBV_4862_4102_CORRHb']=10**(0.4 *klHb *self.hdR['EBV_4862_4102'] )
			self.hdR['EBV_4862_4102_CORRHb_err']= self.hdR['EBV_4862_4102_err'] * n.log(10) * 0.4 * klHb * self.hdR['EBV_4862_4102_CORRHb']
		else :
			self.hdR['BD_4862_4102'] = self.dV
			self.hdR['BD_4862_4102_err'] = self.dV
			self.hdR['EBV_4862_4102'] = self.dV
			self.hdR['EBV_4862_4102_err'] = self.dV
			self.hdR['EBV_4862_4102_CORRO2'] = self.dV
			self.hdR['EBV_4862_4102_CORRO2_err'] = self.dV
			self.hdR['EBV_4862_4102_CORRO3'] = self.dV
			self.hdR['EBV_4862_4102_CORRO3_err'] = self.dV
			self.hdR['EBV_4862_4102_CORRHb'] = self.dV
			self.hdR['EBV_4862_4102_CORRHb_err'] = self.dV

		if self.hdR['H1_4102_flux_nc']>0 and self.hdR['H1_4341_flux_nc']>0 and self.hdR['H1_4102_fluxErr_nc']>0 and self.hdR['H1_4341_fluxErr_nc']>0 :
			self.BDarray[2]=1
			# Balmer decrement : 4341 / 4102
			self.hdR['BD_4102_4341']= self.hdR['H1_4102_flux_nc']/ self.hdR['H1_4341_flux_nc']
			bdc23ErrFrac = ( (self.hdR['H1_4102_fluxErr_nc']/ self.hdR['H1_4102_flux_nc'] )**2 + (self.hdR['H1_4341_fluxErr_nc']/self.hdR['H1_4341_flux_nc'])**2. ) **0.5
			self.hdR['BD_4102_4341_err']=self.hdR['BD_4102_4341'] * bdc23ErrFrac
			# E(B-V) Correction using 4341 / 4102
			self.hdR['EBV_4102_4341'] = -5*n.log10( self.hdR['BD_4102_4341'] * bdc23_ref )/( 2*(5.39 - 5.12))
			self.hdR['EBV_4102_4341_err'] = -5 * bdc23ErrFrac * bdc23_ref /( 2*(5.39 - 5.12)*n.log(10))
			# applied to lines using Calzetti's law
			self.hdR['EBV_4102_4341_CORRO2']=10**(0.4 *self.hdR['EBV_4102_4341'] *klO2)
			self.hdR['EBV_4102_4341_CORRO2_err']= self.hdR['EBV_4102_4341_err'] * n.log(10) * 0.4 * klO2 * self.hdR['EBV_4102_4341_CORRO2']
			self.hdR['EBV_4102_4341_CORRO3'] = 10**(0.4 *self.hdR['EBV_4102_4341'] *klO3)
			self.hdR['EBV_4102_4341_CORRO3_err']= self.hdR['EBV_4102_4341_err'] * n.log(10) * 0.4 * klO3 * self.hdR['EBV_4102_4341_CORRO3']
			self.hdR['EBV_4102_4341_CORRHb']=10**(0.4 *klHb *self.hdR['EBV_4102_4341'] )
			self.hdR['EBV_4102_4341_CORRHb_err']= self.hdR['EBV_4102_4341_err'] * n.log(10) * 0.4 * klHb * self.hdR['EBV_4102_4341_CORRHb']
		else :
			self.hdR['BD_4102_4341'] = self.dV
			self.hdR['BD_4102_4341_err'] = self.dV
			self.hdR['EBV_4102_4341'] = self.dV
			self.hdR['EBV_4102_4341_err'] = self.dV
			self.hdR['EBV_4102_4341_CORRO2'] = self.dV
			self.hdR['EBV_4102_4341_CORRO2_err'] = self.dV
			self.hdR['EBV_4102_4341_CORRO3'] = self.dV
			self.hdR['EBV_4102_4341_CORRO3_err'] = self.dV
			self.hdR['EBV_4102_4341_CORRHb'] = self.dV
			self.hdR['EBV_4102_4341_CORRHb_err'] = self.dV

		# if BD computation succeeded, we can compute instrinsic quantities
		cor_names = n.array(['EBV_4862_4341_CORR','EBV_4862_4102_CORR', 'EBV_4102_4341_CORR'])
		if len((self.BDarray==1).nonzero()[0]>=1):
			name = cor_names[(self.BDarray==1)][0]
			# intrinsic fluxes O3 4960
			self.hdR['flux_O3_4960_intrinsic'] = self.hdR['O3_4960_flux_nc']/ self.hdR[name+'O3']
			self.hdR['flux_O3_4960_intrinsic_err']= self.hdR['flux_O3_4960_intrinsic'] * ((self.hdR['O3_4960_fluxErr_nc']/ self.hdR['O3_4960_flux_nc'] )**2.+ (self.hdR[name+'O3_err']/ self.hdR[name+'O3'] )**2.)**0.5
			# intrinsic fluxes O3 5007
			self.hdR['flux_O3_5007_intrinsic'] = self.hdR['O3_5007_flux_nc']/ self.hdR[name+'O3']
			self.hdR['flux_O3_5007_intrinsic_err']=self.hdR['flux_O3_5007_intrinsic'] * ((self.hdR['O3_5007_fluxErr_nc']/ self.hdR['O3_5007_flux_nc'] )**2.+ (self.hdR[name+'O3_err']/ self.hdR[name+'O3'] )**2.)**0.5
			# intrinsic fluxes O2
			self.hdR['flux_O2_3728_intrinsic']= self.hdR['O2_3728_flux_nc']/ self.hdR[name+'O2']
			self.hdR['flux_O2_3728_intrinsic_err']=self.hdR['flux_O2_3728_intrinsic'] * ((self.hdR['O2_3728_fluxErr_nc']/ self.hdR['O2_3728_flux_nc'] )**2.+ (self.hdR[name+'O2_err']/ self.hdR[name+'O2'] )**2.)**0.5
			# intrinsic fluxes Hb
			self.hdR['flux_H1_4862_intrinsic'] =self.hdR['H1_4862_flux_nc'] / self.hdR[name+'Hb']
			self.hdR['flux_H1_4862_intrinsic_err'] =self.hdR['flux_H1_4862_intrinsic'] * ((self.hdR['H1_4862_fluxErr_nc']/ self.hdR['H1_4862_flux_nc'])**2.+ (self.hdR[name+'Hb_err'] /self.hdR[name+'Hb'] )**2.)**0.5
			# deduce SFR from [OII] 
			self.hdR['SFR_O2_3728'] = 10**(0.27) * 10**(-41) * self.hdR['flux_O2_3728_intrinsic'] * self.sphereCM.value
			self.hdR['SFR_O2_3728_err'] = self.hdR['SFR_O2_3728'] * self.hdR['flux_O2_3728_intrinsic_err'] / self.hdR['flux_O2_3728_intrinsic']
			# deduce SFR from Hbeta
			self.hdR['SFR_H1_4862'] = 10**(0.58) * 10**(-41) * self.hdR['flux_H1_4862_intrinsic'] * self.sphereCM.value
			self.hdR['SFR_H1_4862_err'] = self.hdR['SFR_H1_4862'] * self.hdR['flux_H1_4862_intrinsic_err'] / self.hdR['flux_H1_4862_intrinsic']
		else :
			self.hdR['flux_O3_4960_intrinsic'] = self.dV
			self.hdR['flux_O3_4960_intrinsic_err'] = self.dV
			self.hdR['flux_O3_5007_intrinsic'] = self.dV
			self.hdR['flux_O3_5007_intrinsic_err'] = self.dV
			self.hdR['flux_O2_3728_intrinsic'] = self.dV
			self.hdR['flux_O2_3728_intrinsic_err'] = self.dV
			self.hdR['flux_H1_4862_intrinsic'] = self.dV
			self.hdR['flux_H1_4862_intrinsic_err'] = self.dV
			self.hdR['SFR_O2_3728'] = self.dV
			self.hdR['SFR_O2_3728_err'] = self.dV
			self.hdR['SFR_H1_4862'] = self.dV
			self.hdR['SFR_H1_4862_err'] = self.dV

		# computes O32
		if self.hdR['O3_4960_flux_nc']>0 and self.hdR['O3_5007_flux_nc']>0 and self.hdR['O2_3728_flux_nc']>0 and self.hdR['O3_4960_fluxErr_nc']>0 and self.hdR['O3_5007_fluxErr_nc']>0 and self.hdR['O2_3728_fluxErr_nc'] >0 :
			self.hdR['O32'] = (self.hdR['O3_4960_flux_nc']+ self.hdR['O3_5007_flux_nc'])/ self.hdR['O2_3728_flux_nc']
			O32ErrFrac =  ( ((self.hdR['O3_4960_fluxErr_nc']+ self.hdR['O3_5007_fluxErr_nc'])/ (self.hdR['O3_4960_flux_nc'] +self.hdR['O3_5007_flux_nc']))**2. + (self.hdR['O2_3728_fluxErr_nc'] /self.hdR['O2_3728_flux_nc'] ) **2.)**0.5  
			self.hdR['O32_err'] = self.hdR['O32'] * O32ErrFrac
			# compute R23
			self.hdR['R23'] = (self.hdR['O3_4960_flux_nc']+self.hdR['O3_5007_flux_nc']+ self.hdR['O2_3728_flux_nc'])/self.hdR['H1_4862_flux_nc']
			R23ErrFrac=( ((self.hdR['O3_4960_fluxErr_nc']+ self.hdR['O3_5007_fluxErr_nc']+ self.hdR['O2_3728_fluxErr_nc']) / (self.hdR['O3_4960_flux_nc']+ self.hdR['O3_5007_flux_nc']+ self.hdR['O2_3728_flux_nc']))**2. + (self.hdR['H1_4862_fluxErr_nc']/ self.hdR['H1_4862_flux_nc'])**2.)**0.5  
			self.hdR['R23_err'] = self.hdR['R23'] * R23ErrFrac
			# 12 log(O/H) with Tremonti 04 estimator
			if self.hdR['R23']>0:
				self.hdR['12logOH_tremonti04'] = 9.185-0.313*n.log10(self.hdR['R23']) - 0.264 *n.log10(self.hdR['R23'])**2 - 0.321 *n.log10(self.hdR['R23'])**3
				self.hdR['12logOH_tremonti04_err'] = -0.313* R23ErrFrac / n.log(10) - 0.264 *2 * R23ErrFrac / n.log(10) * n.log10(self.hdR['R23']) - 0.321 * 3* R23ErrFrac / n.log(10) * n.log10(self.hdR['R23'])**2

		else :
			self.hdR['O32'] = self.dV
			self.hdR['O32_err'] = self.dV
			self.hdR['R23'] = self.dV
			self.hdR['R23_err'] = self.dV
			self.hdR['12logOH_tremonti04'] = self.dV
			self.hdR['12logOH_tremonti04_err'] = self.dV

		if self.hdR['flux_O3_4960_intrinsic']>0 and self.hdR['flux_O3_5007_intrinsic']>0 and self.hdR['flux_O2_3728_intrinsic']>0 and self.hdR['flux_H1_4862_intrinsic']>0 and self.hdR['flux_O2_3728_intrinsic_err']>0 and self.hdR['flux_O3_5007_intrinsic_err'] >0 :
			# compute R23 intrinsic
			self.hdR['R23_intrinsic'] = (self.hdR['flux_O3_4960_intrinsic']+ self.hdR['flux_O3_5007_intrinsic']+ self.hdR['flux_O2_3728_intrinsic']) /self.hdR['flux_H1_4862_intrinsic']
			R23ErrFrac_intrinsic=( ((self.hdR['flux_O3_4960_intrinsic_err']+ self.hdR['flux_O3_5007_intrinsic_err'] + self.hdR['flux_O2_3728_intrinsic_err']) / (self.hdR['flux_O3_4960_intrinsic'] + self.hdR['flux_O3_5007_intrinsic']+ self.hdR['flux_O2_3728_intrinsic'])) **2. + (self.hdR['flux_H1_4862_intrinsic_err']/ self.hdR['flux_H1_4862_intrinsic']) **2.)**0.5  
			self.hdR['R23_intrinsic_err'] = self.hdR['R23_intrinsic'] * R23ErrFrac_intrinsic
			if self.hdR['R23_intrinsic']>0:
				self.hdR['12logOH_tremonti04_intrinsic'] = 9.185-0.313* n.log10(self.hdR['R23_intrinsic']) - 0.264 *n.log10(self.hdR['R23_intrinsic'])**2 - 0.321 *n.log10(self.hdR['R23_intrinsic'])**3
				self.hdR['12logOH_tremonti04_intrinsic_err'] = -0.313* R23ErrFrac_intrinsic / n.log(10) - 0.264 *2 * R23ErrFrac_intrinsic / n.log(10) * n.log10(self.hdR['R23_intrinsic']) - 0.321 * 3* R23ErrFrac_intrinsic / n.log(10) * n.log10(self.hdR['R23_intrinsic'])**2

			# 12 log OH with O2(3728)/Hbeta
			OpH = self.hdR['O2_3728_flux_nc']/self.hdR['H1_4862_flux_nc']
			OpHErrFrac = ( (self.hdR['O2_3728_fluxErr_nc']/self.hdR['O2_3728_flux_nc'])**2. + (self.hdR['H1_4862_fluxErr_nc']/self.hdR['H1_4862_flux_nc'])**2.)**0.5  
			OpHErr = OpH * OpHErrFrac
			if OpH>0:
				self.hdR['12logO2H'] = n.log10(OpH) + 7.637
				self.hdR['12logO2H_err'] = OpHErrFrac/n.log(10.)

			# 12 log OH with O3(4960+5007)/Hbeta
			O3H = (self.hdR['O3_4960_flux_nc']+ self.hdR['O3_5007_flux_nc'])/ self.hdR['H1_4862_flux_nc']
			O3HErrFrac =  ( ((self.hdR['O3_4960_fluxErr_nc']+self.hdR['O3_5007_fluxErr_nc'])/ (self.hdR['O3_4960_flux_nc']+ self.hdR['O3_5007_flux_nc']))**2. + (self.hdR['H1_4862_fluxErr_nc']/ self.hdR['H1_4862_flux_nc'])**2.)**0.5  
			O3HErr = O3H* O3HErrFrac
			if O3H>0:
				self.hdR['12logO3H'] = n.log10( O3H )+7.437
				self.hdR['12logO3H_err'] = O3HErrFrac / n.log(10.)

			# 12 log OH with O3(5007)/Hbeta
			O35H = (self.hdR['O3_5007_flux_nc'])/self.hdR['H1_4862_flux_nc']
			O35HErrFrac =  ( (self.hdR['O3_5007_fluxErr_nc']/self.hdR['O3_5007_flux_nc'])**2. + (self.hdR['H1_4862_fluxErr_nc']/self.hdR['H1_4862_flux_nc'])**2.)**0.5  
			O35HErr = O35H* O35HErrFrac
			if O35H>0:
				self.hdR['12logO3_5007_H'] = n.log10( O35H )
				self.hdR['12logO3_5007_H_err'] = O35HErrFrac / n.log(10.)

		else :
			self.hdR['R23_intrinsic'] = self.dV
			self.hdR['R23_intrinsic_err'] = self.dV
			self.hdR['12logOH_tremonti04_intrinsic'] = self.dV
			self.hdR['12logOH_tremonti04_intrinsic_err'] = self.dV
			self.hdR['12logO2H'] = self.dV
			self.hdR['12logO2H_err'] = self.dV
			self.hdR['12logO3H'] = self.dV
			self.hdR['12logO3H_err'] = self.dV
			self.hdR['12logO3_5007_H'] = self.dV
			self.hdR['12logO3_5007_H_err'] = self.dV

	def save_spectrum(self):
		"""
		Saves the stack spectrum, the model and derived quantities in a single fits file with different hdus.
		"""
		wavelength = fits.Column(name="wavelength",format="D", unit="Angstorm", array= 			self.wlLineSpectrum)
		flux = fits.Column(name="flux",format="D", unit="Angstorm", array= 			self.flLineSpectrum)
		fluxErr = fits.Column(name="fluxErr",format="D", unit="Angstorm", array= 			self.flErrLineSpectrum)
		# new columns
		cols = fits.ColDefs([wavelength, flux, fluxErr])
		tbhdu = fits.BinTableHDU.from_columns(cols)
		# previous file
		prihdu = fits.PrimaryHDU(header=self.hdR)
		thdulist = fits.HDUList([prihdu, self.hdu1, self.hdu2, tbhdu])
		outPutFileName = self.stack_model_file[:-5] + "-modeled.fits"
		outFile = np.core.defchararray.replace(outPutFileName, "fits", "model").item()
		os.system('rm '+outFile)
		thdulist.writeto(outFile)

	def plot_fit(self):
		"""
		Plots the fit."""

		age ='age = ' +  str(n.round( 10**self.hdu2.header['light_age'] ,3))+ '+('+ str(n.round( 10**self.hdu2.header['light_age_up']-10**self.hdu2.header['light_age'] ,3)) +')-('+str(n.round( 10**self.hdu2.header['light_age']-10**self.hdu2.header['light_age_low'] ,3))+') Gyr'

		metallicity = 'log(Z/Zsun) = ' + str(n.round( self.hdu2.header['light_metallicity'] ,3))+ '+('+ str(n.round( self.hdu2.header['light_metallicity_up'] - self.hdu2.header['light_metallicity'] ,3)) +')-('+str(n.round( self.hdu2.header['light_metallicity'] - self.hdu2.header['light_metallicity_low'] ,3))+')'

		mass = 'log(M/Msun) = ' + str(n.round( self.hdu2.header['stellar_mass'] ,3))+ '+('+ str(n.round( self.hdu2.header['stellar_mass_up'] - self.hdu2.header['stellar_mass'] ,3)) +')-('+str(n.round( self.hdu2.header['stellar_mass'] - self.hdu2.header['stellar_mass_low'] ,3))+')'

		fig = p.figure(0,(10,10))
		fig.subplots_adjust(hspace=0.5,wspace=0.5)
		# defines the first panel
		fig.add_subplot(4,1,1)
		p.plot(self.stack.x, self.stack.y,'r',label='galaxy',rasterized=True,lw=0.5)
		p.plot(self.wlLineSpectrum, self.model(self.wlLineSpectrum),'b', label='model',rasterized=True)
		p.legend(loc=4)
		#p.xlabel('wavelength')
		p.ylabel(r'$f_\lambda$')
		#p.yscale('log')
		p.xlim((n.min(self.wlLineSpectrum),n.max(self.wlLineSpectrum)))
		p.ylim((0,2*n.max(self.model(self.wlLineSpectrum))))
		p.grid()

		fig.add_subplot(4,1,2)
		p.plot(self.wlLineSpectrum, self.flLineSpectrum, 'k',label='galaxy - model',rasterized=True) #, flErrLineSpectrum
		p.plot(self.wlLineSpectrum, self.stackErr(self.wlLineSpectrum),  'r' , label=' error on data ',rasterized=True)
		p.legend(loc=4)
		#p.xlabel('wavelength')
		p.ylabel(r'$f_\lambda$')
		#p.yscale('log')
		p.xlim((n.min(self.wlLineSpectrum),n.max(self.wlLineSpectrum)))
		p.ylim((-1e-17,1e-17))
		p.grid()

		fig.add_subplot(4,1,3)
		p.plot(self.wlLineSpectrum, self.fl_frac_LineSpectrum, 'k',label='galaxy/model',rasterized=True) #, flErrLineSpectrum
		p.plot(self.wlLineSpectrum, n.ones_like(self.stack(self.wlLineSpectrum)) + self.stackErr(self.wlLineSpectrum)/self.stack(self.wlLineSpectrum),  'r' , label=' frac error on data ',rasterized=True)
		p.legend(loc=4)
		#p.xlabel('wavelength')
		p.ylabel(r'$f_\lambda$')
		#p.yscale('log')
		p.xlim((n.min(self.wlLineSpectrum),n.max(self.wlLineSpectrum)))
		p.ylim((0,5))
		p.grid()

		fig.add_subplot(4,1,4,frame_on=False)
		p.text(0,1,self.stack_file.split('/')[-1]+ ", redshift="+str(n.round(self.redshift,3)) )
		p.text(0,0.8,r' ' + age)
		p.text(0,0.6,metallicity)
		p.text(0,0.4,mass)
		p.tick_params(axis='both', which='both', left='off', right='off', bottom='off', top='off', labelbottom='off', labelleft='off')

		p.savefig(self.stack_file[:-5] + "-fit.pdf")
		p.clf()

	def plot_spectrum(self):
		"""
		Plots the stack spectrum and the model plus derived quantities.
		"""
		# defines the main frame
		fig = p.figure(0,(10,7))
		fig.subplots_adjust(hspace=0.5,wspace=0.5)
		# defines the first panel
		fig.add_subplot(3,1,1)
		p.plot(self.stack.x, self.stack.y,'r',label='stack',lw=2)
		p.plot(self.model.x, self.model.y,'b', label='model')
		p.plot(self.wlLineSpectrum, self.flLineSpectrum, 'k',label='stack - model') #, self.flErrLineSpectrum
		p.legend(loc=2)
		p.xlabel('wavelength')
		p.ylabel(r'$f_\lambda$')
		p.yscale('log')
		p.xlim((n.min(self.model.x),n.max(self.model.x)))
		p.ylim((1e-18,8e-17))
		p.grid()

		fig.add_subplot(3,1,2)
		p.plot(self.wlLineSpectrum, self.flLineSpectrum, 'k',label='stack - model') #, self.flErrLineSpectrum
		p.legend(loc=2)
		p.xlabel('wavelength')
		p.ylabel(r'$f_\lambda$')
		p.yscale('log')
		p.xlim((n.min(self.model.x),n.max(self.model.x)))
		p.ylim((1e-19,1e-17))
		p.grid()

		fig.add_subplot(3,1,3,frame_on=False)
		p.text(0,0,self.stack_file.split('/')[-1])
		tx="BD_4862_4341 = "+str(n.round(self.hdR["BD_4862_4341"],4))+r'$\pm$'+ str(n.round(self.hdR["BD_4862_4341_err"],4))
		p.text(0,0.2,tx)
		tx="EBV_4862_4341 = "+ str(n.round(self.hdR["EBV_4862_4341"],4))+" pm"+ str(n.round(self.hdR["EBV_4862_4341_err"],4))
		p.text(0,0.4,tx)
		tx="EBV_4862_4341_CORRO2 = "+ str(n.round(self.hdR["EBV_4862_4341_CORRO2"],4))+" pm"+ str(n.round(self.hdR["EBV_4862_4341_CORRO2_err"],4))
		p.text(0,0.6,tx)
		tx="BD_4862_4102 = "+ str(n.round(self.hdR["BD_4862_4102"],4))+" pm"+ str(n.round(self.hdR["BD_4862_4102_err"],4))
		p.text(0,0.8,tx)

		p.savefig(self.stack_file[:-5] + "modeled.pdf")
		p.clf()

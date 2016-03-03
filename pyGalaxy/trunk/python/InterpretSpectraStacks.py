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


class InterpretSpectraStacks:
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
			self.stack_spm_file = n.core.defchararray.replace(self.stack_file[:-5], "data", "fits").item() + "-SPM-MILES.fits"
			self.stack_lintFits_file = n.core.defchararray.replace(self.stack_file[:-5], "data", "model").item() + "-SPM-MILES.fits"
		if self.mode=="STELIB":
			self.stack_spm_file = n.core.defchararray.replace(self.stack_file[:-5], "data", "fits").item() + "-SPM-STELIB.fits"
		
		
		
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
		print " loads the data :"
		print self.hdu1.data.dtype
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
		outFile = n.core.defchararray.replace(outPutFileName, "fits", "model").item()
		os.system('rm '+outFile)
		thdulist.writeto(outFile)

		#ADD line model read
		
		#functions :
		
		#combine multiple stack headers in a single table
		
		#summary plot
		
		

		
		
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

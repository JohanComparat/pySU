"""
.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

General purpose:
................

The class ModelSdssSpectra is dedicated to modelling and extracting information from stacks of spectra.

"""
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p
import os 
import astropy.cosmology as co
cosmo=co.Planck13 #co.FlatLambdaCDM(H0=70,Om0=0.3)
import astropy.units as u
import astropy.io.fits as fits
import numpy as n
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.stats import scoreatpercentile

import astropy.io.fits as fits

from lineListVac import *

allLinesList = n.array([ [Ne3,Ne3_3869,"Ne3_3869","left"], [Ne3,Ne3_3968,"Ne3_3968","left"], [O3,O3_4363,"O3_4363","right"], [O3,O3_4960,"O3_4960","left"], [O3,O3_5007,"O3_5007","right"], [N2,N2_6549,"N2_6549","left"], [N2,N2_6585,"N2_6585","right"], [H1,H1_3970,"H1_3970","right"], [H1,H1_4102,"H1_4102","right"], [H1,H1_4341,"H1_4341","right"], [H1,H1_4862,"H1_4862","left"], [H1,H1_6564,"H1_6564","left"]]) 
# other lines that are optional
# , [S2,S2_6718,"S2_6718","left"], [S2,S2_6732,"S2_6732","right"], [Ar3,Ar3_7137,"Ar3_7137","left"], [H1,H1_1216,"H1_1216","right"]

doubletList = n.array([[O2_3727,"O2_3727",O2_3729,"O2_3729",O2_mean]])

# import the fitting routines
import LineFittingLibrary as lineFit

#O2a=3727.092 
#O2b=3729.875 
#O2=(O2a+O2b)/2.
#Hg=4102.892
#Hd=4341.684
#Hb=4862.683
#O3a=4960.295
#O3b=5008.240
#Ha=6564.61

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

class ModelSpectraStacks:
	"""
	This class fits the emission lines on the continuum-subtracted stack.

	:param spec_file: fits file generated with a LF in a luminosity bin.
	:param cosmo: cosmology class from astropy
	:param sdss_min_wavelength: minimum wavelength considered by firefly (default : 1000) 
	:param sdss_max_wavelength: minimum wavelength considered by firefly (default : 7500)
	:param dV: default value that hold the place (default : -9999.99) 
	"""
	def __init__(self, spec_file, mode="MILES", cosmo=cosmo, sdss_min_wavelength= 1000., sdss_max_wavelength=7500., dV=-9999.99, version="stellarpop-m11-salpeter"):
		self.mode = mode
		
		self.spec_file = spec_file
		self.spec_file_base = os.path.basename(spec_file)
		spl = self.spec_file_base.split('.')[0].split('-')
		self.plate = spl[1]
		self.mjd = spl[2]
		self.fiber = spl[3]			
		self.spec_model_file = os.path.join( os.environ['SDSSDR12_DIR'], version, "stellarpop", self.plate, self.spec_file_base[:-5] + "-SPM-MILES.fits")
		outPutDir = os.path.join( os.environ['SDSSDR12_DIR'], version, "model", self.plate)
		if os.path.isdir(outPutDir)==False:
			os.mkdir(outPutDir)
		
		self.outFile = os.path.join( outPutDir, self.spec_file_base[:-5] + ".model")
		
		self.cosmo = cosmo
		self.sdss_max_wavelength	= sdss_max_wavelength
		self.sdss_min_wavelength	= sdss_min_wavelength
		self.dV = dV
		self.side = ''

		hdus = fits.open(self.spec_file)
		self.hdR = hdus[0].header
		self.hdu1 = hdus[1] 
		self.z = hdus[2].data['Z'][0]
		print "Loads the data"
		self.wl = 10**self.hdu1.data['loglam']
		self.fl = self.hdu1.data['flux']
		self.flErr = self.hdu1.data['ivar']**(-0.5)

		self.stack=interp1d(self.wl,self.fl)
		self.stackErr=interp1d(self.wl,self.flErr)
		
		print "loads model"
		hdus = fits.open(self.spec_model_file)
		
		self.hdu2 = hdus[1]
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
		self.stBlue = (self.wl<=self.sdss_min_wavelength)
		# optical
		self.stOpt = (self.wl<self.sdss_max_wavelength)& (self.wl> self.sdss_min_wavelength) 
		# redder than model
		self.stRed = (self.wl>=self.sdss_max_wavelength) 
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
		print "interpolate model"
		self.mdOK =(self.wlModel>n.min(self.wl))&(self.wlModel<n.max(self.wl)) 

		mdBlue=(self.wlModel<=n.min(self.wl)) # bluer part than data
		mdRed=(self.wlModel>=n.max(self.wl)) # redder part than data
		okRed=(self.wlModel>4650)&(self.wlModel<self.sdss_max_wavelength)

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
		print "fits to the line spectrum"
		lfit  =  lineFit.LineFittingLibrary()
		
		#self.subtract_continuum_model()
		data,h=[],[]
		print O2_3727
		dat_mean,mI,hI=lfit.fit_Line_OIIdoublet_position(self.wlLineSpectrum, self.flLineSpectrum, self.flErrLineSpectrum, a0= O2_3727*(1 + self.z) , lineName="O2_3728", p0_sigma=7,model="gaussian",fitWidth = 20.,DLC=10.)
		data.append(dat_mean)
		h.append(hI)

		for li in allLinesList :
			# measure line properties from the mean weighted stack
			print li[2]
			dat_mean,mI,hI=lfit.fit_Line_position_C0noise(self.wlLineSpectrum, self.flLineSpectrum, self.flErrLineSpectrum, li[1]*(1 + self.z), lineName=li[2], continuumSide=li[3], model="gaussian", p0_sigma=7,fitWidth = 15.,DLC=10.)
			data.append(dat_mean)
			h.append(hI)

		heading="".join(h)
		out=n.hstack((data))
		#print "out", out
		out[n.isnan(out)]=n.ones_like(out[n.isnan(out)])*self.dV
		#output = n.array([ out ])
		#print "----------------", output.T[0], output.T[1], output
		colNames = heading.split()
		#print colNames
		col0 = fits.Column(name=colNames[0],format='D', array= n.array([out.T[0]]))
		col1 = fits.Column(name=colNames[1],format='D', array= n.array([out.T[1]]))
		self.lineSpec_cols  = fits.ColDefs([col0, col1])
		#print self.lineSpec_cols
		#print colNames
		for ll in range(2,len(colNames),1):
			#self.hdR["HIERARCH "+colNames[ll]+"_nc"] = out.T[ll]
			self.lineSpec_cols += fits.Column(name=colNames[ll], format='D', array= n.array([out.T[ll]]) )
		
		#print self.lineSpec_cols
		self.lineSpec_tb_hdu = fits.BinTableHDU.from_columns(self.lineSpec_cols)

			
	def fit_lines_to_fullSpectrum(self):
		"""
		Fits the emission lines on the line spectrum.
		"""
		# interpolates the mean spectra.
		print "fits to full spectrum"
		lfit  =  lineFit.LineFittingLibrary()

		data,h=[],[]
		print O2_3727
		dat_mean,mI,hI=lfit.fit_Line_OIIdoublet_position(self.wl, self.fl, self.flErr, a0= O2_3727*(1 + self.z) , lineName="O2_3728", p0_sigma=7,model="gaussian",fitWidth = 20.,DLC=10.)
		print hI, dat_mean
		data.append(dat_mean)
		h.append(hI)

		for li in allLinesList :
			print li[2]
			# measure line properties from the mean weighted stack
			dat_mean,mI,hI=lfit.fit_Line_position_C0noise(self.wl, self.fl, self.flErr, li[1]*(1 + self.z), lineName=li[2], continuumSide=li[3], model="gaussian", p0_sigma=7,fitWidth = 15.,DLC=10.)
			data.append(dat_mean)
			#print li[2], dat_mean
			h.append(hI)

		heading="".join(h)
		out=n.hstack((data))
		out[n.isnan(out)]=n.ones_like(out[n.isnan(out)])*self.dV
		#output = n.array([ out ])
		#print "----------------", output.T[0], output.T[1], output
		colNames = heading.split()
		#print colNames
		col0 = fits.Column(name=colNames[0],format='D', array= n.array([out.T[0]]))
		col1 = fits.Column(name=colNames[1],format='D', array= n.array([out.T[1]]))
		self.fullSpec_cols  = fits.ColDefs([col0, col1])
		#print colNames
		for ll in range(2,len(colNames),1):
			#self.hdR["HIERARCH "+colNames[ll]+"_nc"] = out.T[ll]
			self.fullSpec_cols += fits.Column(name=colNames[ll], format='D', array= n.array([out.T[ll]]) )
		
		self.fullSpec_tb_hdu = fits.BinTableHDU.from_columns(self.fullSpec_cols)


	def save_spectrum(self):
		"""
		Saves the stack spectrum, the model and derived quantities in a single fits file with different hdus.
		"""
		wavelength = fits.Column(name="wavelength",format="D", unit="Angstrom", array= 			self.wlLineSpectrum)
		flux = fits.Column(name="flux",format="D", unit="Angstrom", array= 			self.flLineSpectrum)
		fluxErr = fits.Column(name="fluxErr",format="D", unit="Angstrom", array= 			self.flErrLineSpectrum)
		# new columns
		cols = fits.ColDefs([wavelength, flux, fluxErr])
		lineSptbhdu = fits.BinTableHDU.from_columns(cols)
		# previous file
		prihdu = fits.PrimaryHDU(header=self.hdR)
		thdulist = fits.HDUList([prihdu, self.hdu1, self.hdu2, lineSptbhdu, self.lineSpec_tb_hdu, self.fullSpec_tb_hdu])
		if os.path.isfile(self.outFile):
			os.remove(self.outFile)
		
		thdulist.writeto(self.outFile)

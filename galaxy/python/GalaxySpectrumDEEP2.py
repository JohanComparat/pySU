"""
.. class:: GalaxySpectrumDEEP2

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class GalaxySpectrumDEEP2 is dedicated to handling DEEP2 spectra

"""

dat_2_month = {
  '01': 'jan',
  '02': 'feb',
  '03': 'mar',
  '04': 'apr',
  '05': 'may',
  '06': 'jun',
  '07': 'jul',
  '08': 'aug',
  '09': 'sep',
  '10': 'oct',
  '11': 'nov',
  '12': 'dec'
  }

from os.path import join
import os
import numpy as n
import astropy.io.fits as fits
import glob
from scipy.optimize import curve_fit
from GalaxySurveyDEEP2 import *
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p
from LineFittingLibrary import *
lfl = LineFittingLibrary()
from filterList import *
from lineListAir import *

class GalaxySpectrumDEEP2:
	"""
	Loads the environement proper to the DEEP2 survey.

	Two modes of operation : flux calibration or line fitting
			
	:param catalog_entry: an entry of the deep2 catalog
	:param survey: survey python class
	:param calibration: if the class is loaded with intention of flux calibrating the DEEP2 data.
	:param lineFits: if the class is loaded with intention of fitting line fluxes on the DEEP2 spectra."""
	def __init__(self,catalog_entry, survey=GalaxySurveyDEEP2(redshift_catalog="zcat.deep2.dr4.v4.fits",calibration = False) ,calibration=False,lineFits=True ):
		self.catalog_entry=catalog_entry
		self.mask=str(self.catalog_entry['MASK'])
		self.slit=str(self.catalog_entry['SLIT']).zfill(3)
		self.objno=str(self.catalog_entry['OBJNO'])

		self.database_dir = '/data42s/comparat/firefly/v1_1_0' #os.environ['DATA_DIR']
		self.deep2_dir = join(self.database_dir,"DEEP2")
		self.deep2_catalog_dir = join(self.deep2_dir,"catalogs")
		self.deep2_spectra_dir = join(self.deep2_dir,"raw_data/spectra")
		self.deep2_fc_tc_spectra_dir = join(self.deep2_dir,"raw_data/flux_calibrated_spectra")
		self.survey = survey

		if calibration :
			#self.path_to_spectrum = glob.glob(join(self.deep2_spectra_dir , self.mask ,'*', '*' + self.objno + '*.fits'))
			date_split = self.catalog_entry['DATE'].split('-')
			path_date = "".join(n.array([date_split[0], dat_2_month[date_split[1]] ,date_split[2]]))
			path_folder = join( self.deep2_spectra_dir, self.mask, path_date)
			self.path_to_spectrum = join(path_folder, 'spec1d.'+str(self.catalog_entry['MASK'])+"."+str(self.catalog_entry['SLIT']).zfill(3)+"."+str(self.catalog_entry['OBJNO'])+".fits" )
			self.path_out_folder = join( self.deep2_fc_tc_spectra_dir, self.mask, path_date)
			self.path_to_out_spectrum = join(self.path_out_folder, 'spec1d.'+str(self.catalog_entry['MASK'])+"."+str(self.catalog_entry['SLIT']).zfill(3)+"."+str(self.catalog_entry['OBJNO'])+".fits" )


		if lineFits :
			self.path_to_spectrum = glob.glob(join(self.deep2_spectra_dir , self.mask, '*', '*' + self.objno + '*_fc_tc.dat'))

		#print( "path to spectrum", self.path_to_spectrum, self.mask, self.objno)

		
	def openObservedSpectrum(self):
		"""Loads the observed spectrum in counts."""
		hdS=fits.open(self.path_to_spectrum)
		self.airmass = hdS[1].header['AIRMASS']
		# blue spectrum
		self.dB=hdS[1].data
		# red pectrum
		self.dR=hdS[2].data
		self.chipNO=(hdS[1].header['CHIPNO']-1)%4
		#print hdS[1].header['CHIPNO']-1, self.chipNO
		hdS.close()
		lb=n.hstack((self.dB['LAMBDA'][0],self.dR['LAMBDA'][0]))
		self.lambdSwitch=n.max(self.dB['LAMBDA'][0])
		self.pixSampled=n.arange(2*4096)[(lb>6000)&(lb<10000)]
		self.lambd=lb[(lb>6000)&(lb<10000)]
		self.spec=n.hstack((self.dB['SPEC'][0],self.dR['SPEC'][0]))[(lb>6000)& (lb<10000)]
		self.ivar=n.hstack((self.dB['IVAR'][0],self.dR['IVAR'][0]))[(lb>6000)& (lb<10000)]
		self.specErr=self.ivar**(-0.5)

	def openObservedSpectrumFC(self):
		"""Loads the observed spectrum in counts.
		"""
		self.wavelength,self.fluxl,self.fluxlErr = n.loadtxt(self.path_to_spectrum , unpack=True )

	def correctQE(self):
		"""Corrects from the quantum efficiency of the chips where the spectrum landed."""
		if self.lambd.max()-self.lambd.min() > 3000 or n.mean(self.lambd)<7300 or n.mean(self.lambd)>8300 :
			print( "cannot QE correct" )

		xravg = 8900
		yravg = 150
		correctionavg = self.survey.paramsEndr[0] + self.survey.paramsEndr[1] * self.lambd
		self.xavg = (self.lambd - xravg)/yravg 
		ok1 =  (self.xavg > 0) & ( self.xavg < 1)
		self.cor2avg = correctionavg*self.xavg + 1*(1-self.xavg)
		ok2=(ok1)&(self.cor2avg>1)
		self.cor2avg[(ok2==False)] = n.ones_like(self.cor2avg[(ok2==False)])

		#npixel=len(self.lambd)
		self.left=(self.lambd<=self.lambdSwitch) # n.arange(4096)
		self.right=(self.lambd>self.lambdSwitch) # n.arange(4096,4096*2,1)

		#xx_b=self.lambd[self.left]
		#xx_r=self.lambd[self.right]

		#corr_b = params[num,0] + params[num,1]*self.lambd[self.left] + params[num,2]*self.lambd[self.left]**2
		#corr_r = params[num+4,0] + params[num+4,1]*self.lambd[self.right] + params[num+4,2]*self.lambd[self.right]**2
		corr_b = 1./( self.survey.params.T[self.chipNO][0] + self.survey.params.T[self.chipNO][1] * self.lambd[self.left] + self.survey.params.T[self.chipNO][2]*self.lambd[self.left]**2 )
		corr_r = 1./( self.survey.params.T[self.chipNO+4][0] + self.survey.params.T[self.chipNO+4][1]* self.lambd[self.right] + self.survey.params.T[self.chipNO+4][2] *self.lambd[self.right]**2 )
		# print( corr_b, corr_r, self.cor2avg)
		# print) "spectrum",self.spec)

		self.specTMP=n.zeros_like(self.spec)
		self.specErrTMP=n.zeros_like(self.specErr)
		self.ivarTMP=n.zeros_like(self.ivar)

		self.specTMP[self.left]=self.spec[self.left]*corr_b
		self.specTMP[self.right]=self.spec[self.right]*corr_r* self.cor2avg[self.right]

		self.specErrTMP[self.left]=self.specErr[self.left]*corr_b
		self.specErrTMP[self.right]=self.specErr[self.right]*corr_r* self.cor2avg[self.right]

		self.ivarTMP[self.left]=self.ivar[self.left]/(corr_b*corr_b)
		self.ivarTMP[self.right]=self.ivar[self.right]/(corr_r*corr_r* self.cor2avg[self.right]*self.cor2avg[self.right] )

		self.specTMP=self.specTMP/self.survey.throughput.y[self.pixSampled]
		self.specErrTMP=self.specErrTMP/self.survey.throughput.y[self.pixSampled]
		self.ivarTMP=self.ivarTMP*self.survey.throughput.y[self.pixSampled]**2

	def correct_telluric_abs(self):
		""" Future function to correct the observed spectra from tellurica absorption bands. Not yet implemented. """
		correction = self.survey.telluric_A_band_fct(self.lambd)**(self.airmass**0.55) * self.survey.telluric_B_band_fct(self.lambd)**(self.airmass**0.55)
		self.specTMP=self.specTMP/correction
		self.specErrTMP=self.specErrTMP/correction
		self.ivarTMP=self.ivarTMP*correction**2

	def fluxCal(self):
		"""Performs the flux calibration of the spectrum by converting counts to flux units with an interpolation between the B and the I band borad band photometry."""
		countr = n.sum(self.specTMP*self.survey.Rresponse(self.lambd))/ n.sum(self.survey.Rresponse( self.lambd))
		counti = n.sum(self.specTMP*self.survey.Iresponse(self.lambd))/n.sum( self.survey.Iresponse( self.lambd))
		# (in erg/s/cm^2/Hz)
		fluxr = 10**((self.catalog_entry['MAGR'] + 48.6)/(-2.5)) 
		fluxi = 10**((self.catalog_entry['MAGI'] + 48.6)/(-2.5))
		fpcr = fluxr / countr
		fpci = fluxi / counti
		effr = 6599.0889
		effi = 8135.4026
		x = [effr, effi]
		y = [fpcr, fpci]
		# print( x, y )
		if y[0]>0 and y[1]>0:
			pfits = curve_fit(self.survey.fun,n.log(x),n.log(y),p0=(-0.01,-68))
			fluxn_corr = n.e**( pfits[0][1] + n.log(self.lambd)*pfits[0][0] )
		elif y[0]>0 and y[1]<0:
			fluxn_corr=fpcr
		elif y[0]<0 and y[1]>0:
			fluxn_corr=fpci
		else :
			return "bad"

		self.fluxn = fluxn_corr * self.specTMP
		self.fluxnErr = fluxn_corr * self.specErrTMP
		self.ivar_fluxn=self.ivarTMP/fluxn_corr**2

		self.fluxl=self.fluxn*299792458.0 / (self.lambd**2 * 10**(-10))
		self.fluxlErr=self.fluxnErr *299792458.0 / (self.lambd**2 * 10**(-10))
		self.slit_correction = fluxn_corr
		print(n.median(self.slit_correction))
		return fluxn_corr


	def writeFCspec(self):
		"""Writes the flux-calibrated spectrum"""
		#print(self.path_out_folder)
		if os.path.isdir(self.path_out_folder)==False:
			#print("test")
			os.system('mkdir -p '+self.path_out_folder)
		
		#ff=open(self.path_to_spectrum[:-5]+"_fc_tc.dat",'w')
		#n.savetxt(ff,n.transpose([self.lambd,self.fluxl,self.fluxlErr]))
		#ff.close()

		prihdr = fits.Header()
		prihdr['FILE']          = os.path.basename(self.path_to_out_spectrum)
		prihdr['MASK']          = self.catalog_entry['MASK'] 
		prihdr['OBJNO']         = self.catalog_entry['OBJNO']  
		prihdr['SLIT']          = self.catalog_entry['SLIT']
		prihdr['RA']            = self.catalog_entry['RA']
		prihdr['DEC']           = self.catalog_entry['DEC']
		prihdr['redshift']      = self.catalog_entry['ZBEST']
		prihdr['flux_cal']      = n.log10(n.median(self.slit_correction))
		prihdu = fits.PrimaryHDU(header=prihdr)

		waveCol = fits.Column(name="wavelength",format="D", unit="Angstrom", array= self.lambd)
		dataCol = fits.Column(name="flux",format="D", unit="1e-17erg/s/cm2/Angstrom", array= self.fluxl)
		errorCol = fits.Column(name="flux_error",format="D", unit="1e-17erg/s/cm2/Angstrom", array= self.fluxlErr)
		
		cols = fits.ColDefs([ waveCol, dataCol, errorCol]) 
		tbhdu = fits.BinTableHDU.from_columns(cols)


		complete_hdus = fits.HDUList([prihdu, tbhdu])
		if os.path.isfile(self.path_to_out_spectrum):
			os.remove(self.path_to_out_spectrum)
		complete_hdus.writeto(self.path_to_out_spectrum)
		print("done: ",self.path_to_out_spectrum)


	def openCalibratedSpectrum(self):
		"""Loads the flux calibrated spectrum in f lambda convention.
		"""
		self.wavelength,self.fluxl,self.fluxlErr= n.loadtxt(self.path_to_spectrum[0], unpack=True)


	def plotFit(self, outputFigureNameRoot, ymin = 1e-19, ymax = 1e-17):
		"""
		Plots spectrum and the line fits in a few figures
		:param outputFigureNameRoot:  path + name to save the plots
		"""
		
		ok = (self.fluxl >0 ) & (self.fluxl > 1.5* self.fluxlErr)
		p.figure(1,(12,4))
		p.axes([0.1,0.2,0.85,0.75])
		p.errorbar(self.wavelength[ok][::4],self.fluxl[ok][::4],yerr = self.fluxlErr[ok][::4], linewidth=1, alpha= 0.4, label='spectrum')
		p.xlabel('wavelength [A]')
		p.ylabel(r'f$_\lambda$ [erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
		p.yscale('log')
		p.xlim((6400,9200))
		p.ylim((ymin, ymax))
		gl = p.legend(loc=0,fontsize=12)
		gl.set_frame_on(False)
		p.savefig( outputFigureNameRoot + "-all.png" )
		p.clf()

		a0_1 = self.catalog_entry['O2_3728_a0a']
		a0_2 = self.catalog_entry['O2_3728_a0b']
		continu= self.catalog_entry['O2_3728_continu']
		aas =n.arange(self.catalog_entry['O2_3728_a0a']-25, self.catalog_entry['O2_3728_a0b']+25,0.05)
		flMod=lambda aa,sigma,F0,sh :continu+ lfl.gaussianLineNC(aa,sigma,(1-sh)*F0,a0_1)+lfl.gaussianLineNC(aa,sigma,sh*F0,a0_2)
		model = flMod(aas, self.catalog_entry['O2_3728_sigma'], self.catalog_entry['O2_3728_flux'], self.catalog_entry['O2_3728_share'] )		
		p.figure(2,(4,4))
		p.axes([0.21,0.2,0.78,0.7])
		p.errorbar(self.wavelength,self.fluxl,yerr = self.fluxlErr)
		p.plot(aas, model,'g',label='model', lw=2)
		p.xlabel('wavelength [A]')
		p.ylabel(r'f$_\lambda$ [erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
		p.yscale('log')
		p.ylim((ymin, ymax))
		p.xlim(( self.catalog_entry['O2_3728_a0a']-25, self.catalog_entry['O2_3728_a0b']+25))
		gl = p.legend(loc=0,fontsize=12)
		gl.set_frame_on(False)
		p.title('[OII] 3727')
		p.savefig( outputFigureNameRoot + "-O2_3728.png")
		p.clf()

		a0 = self.catalog_entry['O3_5007_a0']
		continu= self.catalog_entry['O3_5007_continu']
		aas =n.arange(self.catalog_entry['O3_5007_a0']-25, self.catalog_entry['O3_5007_a0']+25,0.05)
		flMod=lambda aa,sigma,F0: lfl.gaussianLine(aa,sigma,F0,a0,continu)
		model = flMod(aas, self.catalog_entry['O3_5007_sigma'], self.catalog_entry['O3_5007_flux'])
		
		p.figure(2,(4,4))
		p.axes([0.21,0.2,0.78,0.7])
		p.errorbar(self.wavelength,self.fluxl,yerr = self.fluxlErr)
		p.plot(aas, model,'g',label='model', lw =2)
		p.xlabel('wavelength [A]')
		p.ylabel(r'f$_\lambda$ [erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
		p.yscale('log')
		p.ylim((ymin, ymax))
		p.xlim(( self.catalog_entry['O3_5007_a0']-25, self.catalog_entry['O3_5007_a0']+25))
		gl = p.legend(loc=0,fontsize=12)
		gl.set_frame_on(False)
		p.title('[OIII] 5007')
		p.savefig( outputFigureNameRoot + "-O3_5007.png")
		p.clf()

		a0 = self.catalog_entry['H1_4862_a0']
		continu= self.catalog_entry['H1_4862_continu']
		aas =n.arange(self.catalog_entry['H1_4862_a0']-25, self.catalog_entry['H1_4862_a0']+25,0.05)
		flMod=lambda aa,sigma,F0: lfl.gaussianLine(aa,sigma,F0,a0,continu)
		model = flMod(aas, self.catalog_entry['H1_4862_sigma'], self.catalog_entry['H1_4862_flux'])
		
		p.figure(2,(4,4))
		p.axes([0.21,0.2,0.78,0.7])
		p.errorbar(self.wavelength,self.fluxl,yerr = self.fluxlErr)
		p.plot(aas, model,'g',label='model', lw =2)
		p.xlabel('wavelength [A]')
		p.ylabel(r'f$_\lambda$ [erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
		p.yscale('log')
		p.ylim((ymin, ymax))
		p.xlim(( self.catalog_entry['H1_4862_a0']-25, self.catalog_entry['H1_4862_a0']+25))
		gl = p.legend(loc=0,fontsize=12)
		gl.set_frame_on(False)
		p.title(r'H$\beta$')
		p.savefig( outputFigureNameRoot + "-H1_4862.png")
		p.clf()



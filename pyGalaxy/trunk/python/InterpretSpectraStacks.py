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

class InterpretSpectraStacks:
	"""
	This class helps interpret the combination of the fits of the SPM and line fits.

	:param stack_file: fits file generated with a LF in a luminosity bin.
	:param cosmo: cosmology class from astropy
	:param dV: default value that hold the place (default : -9999.99) 
	"""
	def __init__(self, stack_file, mode="MILES", cosmo=cosmo, dV=-9999.99):
		self.stack_file = stack_file
		self.mode = mode
		if self.mode=="MILES":
			self.stack_spm_file = n.core.defchararray.replace(self.stack_file[:-5], "data", "fits").item() + "-SPM-MILES.fits"
			self.stack_lineFits_file = n.core.defchararray.replace(self.stack_file[:-5], "data", "model").item() + "-SPM-MILES-modeled.model"
		if self.mode=="STELIB":
			self.stack_spm_file = n.core.defchararray.replace(self.stack_file[:-5], "data", "fits").item() + "-SPM-STELIB.fits"
			self.stack_lineFits_file = n.core.defchararray.replace(self.stack_file[:-5], "data", "model").item() + "-SPM-STELIB-modeled.model"
		
		self.cosmo = cosmo
		self.dV = dV
		# define global quantities 
		self.redshift = float(self.stack_file.split('-')[2].split('_')[0][1:])
		self.lineWave = int(self.stack_file.split('/')[-1].split('_')[1][:4])
		self.N_in_stack = int(self.stack_file.split('/')[-1].split('_')[4])
		self.R_stack = int(self.stack_file.split('/')[-1].split('_')[6])
		if  self.stack_file.split('/')[-1].split('_')[1].split('-')[1] == "DEEP2R24.2":
			self.survey = int(1)
		if  self.stack_file.split('/')[-1].split('_')[1].split('-')[1] == "VVDSDEEPI24":
			self.survey = int(2)
			
		# opens the stack
		print " loads the stack :"
		self.hduStack = fits.open(self.stack_file)
		# opens the spm model
		print " loads the spm model :"
		self.hduSPM = fits.open(self.stack_spm_file)
		# opens the line model
		print " loads the line model :"
		self.hduLine = fits.open(self.stack_lineFits_file)
		
	def get_table_entry(self):
		header = lineWavelength Survey Redshift L_MIN L_MAX L_MEAN N_in_stack R_stack spm_light_age spm_light_age_err_plus spm_light_age_err_minus spm_light_metallicity spm_light_metallicity_err_plus spm_light_metallicity_err_minus spm_stellar_mass spm_stellar_mass_err_plus spm_stellar_mass_err_minus spm_EBV gp_EBV_4862_4341 gp_EBV_4862_4341_err gp_EBV_4862_4102 gp_EBV_4862_4102_err gp_BD_4102_4341 gp_BD_4102_4341_err gp_SFR_O2_3728 gp_SFR_O2_3728_err gp_SFR_H1_4862 gp_SFR_H1_4862_err gp_12logOH_tremonti04 gp_12logOH_tremonti04_err gp_12logOH_tremonti04_intrinsic gp_12logOH_tremonti04_intrinsic_err 
		table_entry = n.array([self.lineWave, self.survey, self.redshift, self.hduLine[0].header['L_MIN'], self.hduLine[0].header['L_MAX'], self.hduLine[0].header['L_MEAN'], self.N_in_stack, self.R_stack, 10**self.hduSPM[1].header['light_age'], 10**self.hduSPM[1].header['light_age_up']-10**self.hduSPM[1].header['light_age'], 10**self.hduSPM[1].header['light_age']-10**self.hduSPM[1].header['light_age_low'],self.hduSPM[1].header['light_metallicity'], self.hduSPM[1].header['light_metallicity_up'] - self.hduSPM[1].header['light_metallicity'], self.hduSPM[1].header['light_metallicity'] - self.hduSPM[1].header['light_metallicity_low'], self.hduSPM[1].header['stellar_mass'], self.hduSPM[1].header['stellar_mass_up'] - self.hduSPM[1].header['stellar_mass'], self.hduSPM[1].header['stellar_mass'] - self.hduSPM[1].header['stellar_mass_low'], self.hduSPM[1].header['EBV'], self.hduLine[0].header['EBV_4862_4341'], self.hduLine[0].header['EBV_4862_4341_err'], self.hduLine[0].header['EBV_4862_4102'], self.hduLine[0].header['EBV_4862_4102_err'], self.hduLine[0].header['BD_4102_4341'], self.hduLine[0].header['BD_4102_4341_err'], self.hduLine[0].header['SFR_O2_3728'], self.hduLine[0].header['SFR_O2_3728_err'], self.hduLine[0].header['SFR_H1_4862'], self.hduLine[0].header['SFR_H1_4862_err'], self.hduLine[0].header['12logOH_tremonti04'], self.hduLine[0].header['12logOH_tremonti04_err'], self.hduLine[0].header['12logOH_tremonti04_intrinsic'], self.hduLine[0].header['12logOH_tremonti04_intrinsic_err'] ])	
		return table_entry, header
		
	def print_spm_result(self):
		age ='age = ' +  str(n.round( 10**self.hdu2.header['light_age'] ,3))+ '+('+ str(n.round( 10**self.hdu2.header['light_age_up']-10**self.hdu2.header['light_age'] ,3)) +')-('+str(n.round( 10**self.hdu2.header['light_age']-10**self.hdu2.header['light_age_low'] ,3))+') Gyr'
		metallicity = 'log(Z/Zsun) = ' + str(n.round( self.hdu2.header['light_metallicity'] ,3))+ '+('+ str(n.round( self.hdu2.header['light_metallicity_up'] - self.hdu2.header['light_metallicity'] ,3)) +')-('+str(n.round( self.hdu2.header['light_metallicity'] - self.hdu2.header['light_metallicity_low'] ,3))+')'
		mass = 'log(M/Msun) = ' + str(n.round( self.hdu2.header['stellar_mass'] ,3))+ '+('+ str(n.round( self.hdu2.header['stellar_mass_up'] - self.hdu2.header['stellar_mass'] ,3)) +')-('+str(n.round( self.hdu2.header['stellar_mass'] - self.hdu2.header['stellar_mass_low'] ,3))+')'
		return age, metallicity, mass
		

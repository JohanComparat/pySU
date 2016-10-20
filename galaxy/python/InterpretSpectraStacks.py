"""
.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

General purpose:
................

The class InterpretSpectraStacks is dedicated to modelling and extracting information from stacks of spectra.

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


from lineListVac import *
allLinesList = n.array([ [Ne3,Ne3_3869,"Ne3_3869","left"], [O3,O3_4363,"O3_4363","right"], [O3,O3_4960,"O3_4960","left"], [O3,O3_5007,"O3_5007","right"], [N2,N2_6549,"N2_6549","left"], [N2,N2_6585,"N2_6585","right"], [H1,H1_3970,"H1_3970","right"], [H1,H1_4102,"H1_4102","right"], [H1,H1_4341,"H1_4341","right"], [H1,H1_4862,"H1_4862","left"], [H1,H1_6564,"H1_6564","left"]]) 
# other lines that are optional
# , [S2,S2_6718,"S2_6718","left"], [S2,S2_6732,"S2_6732","right"], [Ar3,Ar3_7137,"Ar3_7137","left"], [H1,H1_1216,"H1_1216","right"]

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
	This class helps interpret the combination of the fits of the SPM and line fits.

	:param stack_file: fits file generated with a LF in a luminosity bin.
	:param cosmo: cosmology class from astropy
	:param dV: default value that hold the place (default : -9999.99) 
	"""
	def __init__(self, stack_file, cosmo=cosmo, dV=-9999.99):
		self.stack_file = stack_file
		self.cosmo = cosmo
		self.dV = dV
		# define global quantities 
		self.baseN = os.path.basename(self.stack_file)
		self.redshift = float(self.stack_file.split('-')[2].split('_')[0][1:])
		self.lineWave = int(self.baseN.split('_')[1][:4])
		self.N_in_stack = int(self.baseN.split('_')[4])
		self.R_stack = int(self.baseN.split('_')[6])
		if  self.baseN.split('_')[1].split('-')[1] == "DEEP2R24.2":
			self.survey = int(1)
		if  self.baseN.split('_')[1].split('-')[1] == "VVDSDEEPI24":
			self.survey = int(2)
		
		# opens the stack
		print " loads the stack :", self.stack_file
		hduList = fits.open(self.stack_file)
		self.head = hduList[0]
		self.hduStack = hduList[1]
		wl= self.hduStack.data['wavelength'][(self.hduStack.data['NspectraPerPixel']>0.8 * self.N_in_stack)]
		self.wlmin = n.min(wl)
		self.wlmax = n.max(wl)
		# opens the spm model
		print " loads the spm model :"
		self.hduSPM = hduList[2]
		# opens the line model
		print " loads the line model lineSpec :"
		self.hduLine = hduList[4]
		print " loads the line model fullSpec :"
		self.hduFull = hduList[5]
		
	def get_table_entry_full(self):
		headerA =" lineWavelength Survey Redshift L_MIN L_MAX L_MEAN N_in_stack wl_min wl_max R_stack age_universe age_lightW_mean age_lightW_err_plus age_lightW_err_minus metallicity_lightW_mean metallicity_lightW_mean_err_plus metallicity_lightW_mean_err_minus age_massW_mean age_massW_err_plus age_massW_err_minus metallicity_massW_mean metallicity_massW_mean_err_plus metallicity_massW_mean_err_minus stellar_mass stellar_mass_err_plus stellar_mass_err_minus spm_EBV nComponentsSSP"
		
		table_entry = [self.lineWave, self.survey, self.redshift, self.head.header['L_MIN'], self.head.header['L_MAX'], self.head.header['L_MEAN'], self.N_in_stack, self.wlmin, self.wlmax, self.R_stack, 10**self.hduSPM.header['age_universe'], 10**self.hduSPM.header['age_lightW_mean'], 10**self.hduSPM.header['age_lightW_mean_up']-10**self.hduSPM.header['age_lightW_mean'], 10**self.hduSPM.header['age_lightW_mean']-10**self.hduSPM.header['age_lightW_mean_low'], self.hduSPM.header['metallicity_lightW_mean'], self.hduSPM.header['metallicity_lightW_mean_up'] - self.hduSPM.header['metallicity_lightW_mean'], self.hduSPM.header['metallicity_lightW_mean'] - self.hduSPM.header['metallicity_lightW_mean_low'], 10**self.hduSPM.header['age_massW_mean'], 10**self.hduSPM.header['age_massW_mean_up']-10**self.hduSPM.header['age_massW_mean'], 10**self.hduSPM.header['age_massW_mean']-10**self.hduSPM.header['age_massW_mean_low'], self.hduSPM.header['metallicity_massW_mean'], self.hduSPM.header['metallicity_massW_mean_up'] - self.hduSPM.header['metallicity_massW_mean'], self.hduSPM.header['metallicity_massW_mean'] - self.hduSPM.header['metallicity_massW_mean_low'], self.hduSPM.header['stellar_mass_mean'], self.hduSPM.header['stellar_mass_mean_up'] - self.hduSPM.header['stellar_mass_mean'], self.hduSPM.header['stellar_mass_mean'] - self.hduSPM.header['stellar_mass_mean_low'], self.hduSPM.header['EBV'], self.hduSPM.header['ssp_number']]
		#print self.hduSPM.header
		for iii in n.arange(self.hduSPM.header['ssp_number']):
			table_entry.append( self.hduSPM.header['stellar_mass_ssp_'+str(iii)] )
			table_entry.append( self.hduSPM.header['age_ssp_'+str(iii)] )
			table_entry.append( self.hduSPM.header['metal_ssp_'+str(iii)] )
			table_entry.append( self.hduSPM.header['SFR_ssp_'+str(iii)] )
			table_entry.append( self.hduSPM.header['weightMass_ssp_'+str(iii)] )
			table_entry.append( self.hduSPM.header['weightLight_ssp_'+str(iii)] )
			headerA += ' stellar_mass_ssp_'+str(iii) + ' age_ssp_'+str(iii) + ' metal_ssp_'+str(iii) + ' SFR_ssp_'+str(iii) + ' weightMass_ssp_'+str(iii) + ' weightLight_ssp_'+str(iii)
		
		if self.hduSPM.header['ssp_number']<8 :
			for iii in n.arange(self.hduSPM.header['ssp_number'], 8, 1):
				table_entry.append([0., 0., 0., 0., 0., 0.])
				headerA += ' stellar_mass_ssp_'+str(iii) + ' age_ssp_'+str(iii) + ' metal_ssp_'+str(iii) + ' SFR_ssp_'+str(iii) + ' weightMass_ssp_'+str(iii) + ' weightLight_ssp_'+str(iii)

		table_entry = n.array( n.hstack((table_entry)) )
		#print table_entry.shape
		return n.hstack((table_entry, self.hduFull.data[0])), headerA
		
	def get_table_entry_line(self):
		headerA =" lineWavelength Survey Redshift L_MIN L_MAX L_MEAN N_in_stack wl_min wl_max R_stack age_universe age_lightW_mean age_lightW_err_plus age_lightW_err_minus metallicity_lightW_mean metallicity_lightW_mean_err_plus metallicity_lightW_mean_err_minus age_massW_mean age_massW_err_plus age_massW_err_minus metallicity_massW_mean metallicity_massW_mean_err_plus metallicity_massW_mean_err_minus stellar_mass stellar_mass_err_plus stellar_mass_err_minus spm_EBV nComponentsSSP"
		
		table_entry = [self.lineWave, self.survey, self.redshift, self.head.header['L_MIN'], self.head.header['L_MAX'], self.head.header['L_MEAN'], self.N_in_stack, self.wlmin, self.wlmax, self.R_stack, 10**self.hduSPM.header['age_universe'], 10**self.hduSPM.header['age_lightW_mean'], 10**self.hduSPM.header['age_lightW_mean_up']-10**self.hduSPM.header['age_lightW_mean'], 10**self.hduSPM.header['age_lightW_mean']-10**self.hduSPM.header['age_lightW_mean_low'], self.hduSPM.header['metallicity_lightW_mean'], self.hduSPM.header['metallicity_lightW_mean_up'] - self.hduSPM.header['metallicity_lightW_mean'], self.hduSPM.header['metallicity_lightW_mean'] - self.hduSPM.header['metallicity_lightW_mean_low'], 10**self.hduSPM.header['age_massW_mean'], 10**self.hduSPM.header['age_massW_mean_up']-10**self.hduSPM.header['age_massW_mean'], 10**self.hduSPM.header['age_massW_mean']-10**self.hduSPM.header['age_massW_mean_low'], self.hduSPM.header['metallicity_massW_mean'], self.hduSPM.header['metallicity_massW_mean_up'] - self.hduSPM.header['metallicity_massW_mean'], self.hduSPM.header['metallicity_massW_mean'] - self.hduSPM.header['metallicity_massW_mean_low'], self.hduSPM.header['stellar_mass_mean'], self.hduSPM.header['stellar_mass_mean_up'] - self.hduSPM.header['stellar_mass_mean'], self.hduSPM.header['stellar_mass_mean'] - self.hduSPM.header['stellar_mass_mean_low'], self.hduSPM.header['EBV'], self.hduSPM.header['ssp_number']]
		#print self.hduSPM.header
		for iii in n.arange(self.hduSPM.header['ssp_number']):
			table_entry.append( self.hduSPM.header['stellar_mass_ssp_'+str(iii)] )
			table_entry.append( self.hduSPM.header['age_ssp_'+str(iii)] )
			table_entry.append( self.hduSPM.header['metal_ssp_'+str(iii)] )
			table_entry.append( self.hduSPM.header['SFR_ssp_'+str(iii)] )
			table_entry.append( self.hduSPM.header['weightMass_ssp_'+str(iii)] )
			table_entry.append( self.hduSPM.header['weightLight_ssp_'+str(iii)] )
			headerA += ' stellar_mass_ssp_'+str(iii) + ' age_ssp_'+str(iii) + ' metal_ssp_'+str(iii) + ' SFR_ssp_'+str(iii) + ' weightMass_ssp_'+str(iii) + ' weightLight_ssp_'+str(iii)
		
		if self.hduSPM.header['ssp_number']<8 :
			for iii in n.arange(self.hduSPM.header['ssp_number'], 8, 1):
				table_entry.append([0., 0., 0., 0., 0., 0.])
				headerA += ' stellar_mass_ssp_'+str(iii) + ' age_ssp_'+str(iii) + ' metal_ssp_'+str(iii) + ' SFR_ssp_'+str(iii) + ' weightMass_ssp_'+str(iii) + ' weightLight_ssp_'+str(iii)

		table_entry = n.array( n.hstack((table_entry)) )
		#print table_entry.shape
		# headerB = "gp_EBV_4862_4341 gp_EBV_4862_4341_err gp_EBV_4862_4102 gp_EBV_4862_4102_err gp_BD_4102_4341 gp_BD_4102_4341_err gp_SFR_O2_3728 gp_SFR_O2_3728_err gp_SFR_H1_4862 gp_SFR_H1_4862_err gp_12logOH_tremonti04 gp_12logOH_tremonti04_err gp_12logOH_tremonti04_intrinsic gp_12logOH_tremonti04_intrinsic_err "
		# table_entry = n.array([ self.hduLine.header['EBV_4862_4341'], self.hduLine.header['EBV_4862_4341_err'], self.hduLine.header['EBV_4862_4102'], self.hduLine.header['EBV_4862_4102_err'], self.hduLine.header['BD_4102_4341'], self.hduLine.header['BD_4102_4341_err'], self.hduLine.header['SFR_O2_3728'], self.hduLine.header['SFR_O2_3728_err'], self.hduLine.header['SFR_H1_4862'], self.hduLine.header['SFR_H1_4862_err'], self.hduLine.header['12logOH_tremonti04'], self.hduLine.header['12logOH_tremonti04_err'], self.hduLine.header['12logOH_tremonti04_intrinsic'], self.hduLine.header['12logOH_tremonti04_intrinsic_err'] ])	
		return n.hstack((table_entry, self.hduLine.data[0])), headerA
		
	def print_spm_result(self):
		age ='age = ' +  str(n.round( 10**self.hdu2.header['light_age'] ,3))+ '+('+ str(n.round( 10**self.hdu2.header['age_lightW_mean_up']-10**self.hdu2.header['light_age'] ,3)) +')-('+str(n.round( 10**self.hdu2.header['age_lightW_mean']-10**self.hdu2.header['age_lightW_mean_low'] ,3))+') Gyr'
		metallicity = 'log(Z/Zsun) = ' + str(n.round( self.hdu2.header['metallicity_lightW_mean'] ,3))+ '+('+ str(n.round( self.hdu2.header['metallicity_lightW_mean_up'] - self.hdu2.header['metallicity_lightW_mean'] ,3)) +')-('+str(n.round( self.hdu2.header['metallicity_lightW_mean'] - self.hdu2.header['metallicity_lightW_mean_low'] ,3))+')'
		mass = 'log(M/Msun) = ' + str(n.round( self.hdu2.header['stellar_mass_mean'] ,3))+ '+('+ str(n.round( self.hdu2.header['stellar_mass_mean_up'] - self.hdu2.header['stellar_mass_mean'] ,3)) +')-('+str(n.round( self.hdu2.header['stellar_mass_mean'] - self.hdu2.header['stellar_mass_mean_low'] ,3))+')'
		return age, metallicity, mass
		
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
			self.hdR['HIERARCH BD_4862_4341']=self.hdR['H1_4341_flux_nc']/ self.hdR['H1_4862_flux_nc']
			bdc1ErrFrac = ( (self.hdR['H1_4862_fluxErr_nc']/ self.hdR['H1_4862_flux_nc'])**2 + (self.hdR['H1_4341_fluxErr_nc']/ self.hdR['H1_4341_flux_nc'])**2. ) **0.5
			self.hdR['HIERARCH BD_4862_4341_err']= self.hdR['BD_4862_4341'] * bdc1ErrFrac
			# E(B-V) Correction using 4862 / 4341
			self.hdR['HIERARCH EBV_4862_4341'] = -5*n.log10(self.hdR['BD_4862_4341'] * bdc1_ref) / (2* (5.12 - 4.6))
			self.hdR['HIERARCH EBV_4862_4341_err']= -5 * bdc1ErrFrac * bdc1_ref/(2*(5.12-4.6)*n.log(10))
			# applied to emission lines using Calzetti's law
			self.hdR['HIERARCH EBV_4862_4341_CORRO2']=10**(0.4 * self.hdR['EBV_4862_4341'] *klO2)
			self.hdR['HIERARCH EBV_4862_4341_CORRO2_err']= self.hdR['EBV_4862_4341_err'] * n.log(10) * 0.4 * klO2 * self.hdR['EBV_4862_4341_CORRO2']
			self.hdR['HIERARCH EBV_4862_4341_CORRO3']= 10**(0.4 *self.hdR['EBV_4862_4341'] *klO3)
			self.hdR['HIERARCH EBV_4862_4341_CORRO3_err']= self.hdR['EBV_4862_4341_err'] * n.log(10) * 0.4 * klO3 * self.hdR['EBV_4862_4341_CORRO3']
			self.hdR['HIERARCH EBV_4862_4341_CORRHb']=10**(0.4 *klHb *self.hdR['EBV_4862_4341'] )
			self.hdR['HIERARCH EBV_4862_4341_CORRHb_err']= self.hdR['EBV_4862_4341_err'] * n.log(10) * 0.4 * klHb * self.hdR['EBV_4862_4341_CORRHb']
		else :
			self.hdR['HIERARCH BD_4862_4341'] = self.dV
			self.hdR['HIERARCH BD_4862_4341_err'] = self.dV
			self.hdR['HIERARCH EBV_4862_4341'] = self.dV
			self.hdR['HIERARCH EBV_4862_4341_err'] = self.dV
			self.hdR['HIERARCH EBV_4862_4341_CORRO2'] = self.dV
			self.hdR['HIERARCH EBV_4862_4341_CORRO2_err'] = self.dV
			self.hdR['HIERARCH EBV_4862_4341_CORRO3'] = self.dV
			self.hdR['HIERARCH EBV_4862_4341_CORRO3_err'] = self.dV
			self.hdR['HIERARCH EBV_4862_4341_CORRHb'] = self.dV
			self.hdR['HIERARCH EBV_4862_4341_CORRHb_err'] = self.dV

		if self.hdR['H1_4102_flux_nc']>0 and self.hdR['H1_4862_flux_nc']>0 and self.hdR['H1_4102_fluxErr_nc']>0 and self.hdR['H1_4862_fluxErr_nc']>0 :
			self.BDarray[1]=1
			# Balmer decrement : 4862 / 4102
			self.hdR['HIERARCH BD_4862_4102'] = self.hdR['H1_4102_flux_nc']/ self.hdR['H1_4862_flux_nc']
			bdc2ErrFrac = ( (self.hdR['H1_4862_fluxErr_nc']/ self.hdR['H1_4862_flux_nc'] )**2 + (self.hdR['H1_4102_fluxErr_nc']/self.hdR['H1_4102_flux_nc'])**2. ) **0.5
			self.hdR['HIERARCH BD_4862_4102_err'] = self.hdR['BD_4862_4102']* bdc2ErrFrac
			# E(B-V) Correction using 4862 / 4341
			self.hdR['HIERARCH EBV_4862_4102'] = -5*n.log10( self.hdR['BD_4862_4102'] * bdc2_ref )/( 2*(5.39-4.6))
			self.hdR['HIERARCH EBV_4862_4102_err'] = -5 * bdc2ErrFrac * bdc2_ref /(2*( 5.39 - 4.6)*n.log(10))
			# applied to emission lines using Calzetti's law
			self.hdR['HIERARCH EBV_4862_4102_CORRO2']=10**(0.4 *self.hdR['EBV_4862_4102'] *klO2)
			self.hdR['HIERARCH EBV_4862_4102_CORRO2_err']= self.hdR['EBV_4862_4102_err'] * n.log(10) * 0.4 * klO2 * self.hdR['EBV_4862_4102_CORRO2']
			self.hdR['HIERARCH EBV_4862_4102_CORRO3']= 10**(0.4 *self.hdR['EBV_4862_4102'] *klO3)
			self.hdR['HIERARCH EBV_4862_4102_CORRO3_err']= self.hdR['EBV_4862_4102_err'] * n.log(10) * 0.4 * klO3 * self.hdR['EBV_4862_4102_CORRO3']
			self.hdR['HIERARCH EBV_4862_4102_CORRHb']=10**(0.4 *klHb *self.hdR['EBV_4862_4102'] )
			self.hdR['HIERARCH EBV_4862_4102_CORRHb_err']= self.hdR['EBV_4862_4102_err'] * n.log(10) * 0.4 * klHb * self.hdR['EBV_4862_4102_CORRHb']
		else :
			self.hdR['HIERARCH BD_4862_4102'] = self.dV
			self.hdR['HIERARCH BD_4862_4102_err'] = self.dV
			self.hdR['HIERARCH EBV_4862_4102'] = self.dV
			self.hdR['HIERARCH EBV_4862_4102_err'] = self.dV
			self.hdR['HIERARCH EBV_4862_4102_CORRO2'] = self.dV
			self.hdR['HIERARCH EBV_4862_4102_CORRO2_err'] = self.dV
			self.hdR['HIERARCH EBV_4862_4102_CORRO3'] = self.dV
			self.hdR['HIERARCH EBV_4862_4102_CORRO3_err'] = self.dV
			self.hdR['HIERARCH EBV_4862_4102_CORRHb'] = self.dV
			self.hdR['HIERARCH EBV_4862_4102_CORRHb_err'] = self.dV

		if self.hdR['H1_4102_flux_nc']>0 and self.hdR['H1_4341_flux_nc']>0 and self.hdR['H1_4102_fluxErr_nc']>0 and self.hdR['H1_4341_fluxErr_nc']>0 :
			self.BDarray[2]=1
			# Balmer decrement : 4341 / 4102
			self.hdR['HIERARCH BD_4102_4341']= self.hdR['H1_4102_flux_nc']/ self.hdR['H1_4341_flux_nc']
			bdc23ErrFrac = ( (self.hdR['H1_4102_fluxErr_nc']/ self.hdR['H1_4102_flux_nc'] )**2 + (self.hdR['H1_4341_fluxErr_nc']/self.hdR['H1_4341_flux_nc'])**2. ) **0.5
			self.hdR['HIERARCH BD_4102_4341_err']=self.hdR['BD_4102_4341'] * bdc23ErrFrac
			# E(B-V) Correction using 4341 / 4102
			self.hdR['HIERARCH EBV_4102_4341'] = -5*n.log10( self.hdR['BD_4102_4341'] * bdc23_ref )/( 2*(5.39 - 5.12))
			self.hdR['HIERARCH EBV_4102_4341_err'] = -5 * bdc23ErrFrac * bdc23_ref /( 2*(5.39 - 5.12)*n.log(10))
			# applied to lines using Calzetti's law
			self.hdR['HIERARCH EBV_4102_4341_CORRO2']=10**(0.4 *self.hdR['EBV_4102_4341'] *klO2)
			self.hdR['HIERARCH EBV_4102_4341_CORRO2_err']= self.hdR['EBV_4102_4341_err'] * n.log(10) * 0.4 * klO2 * self.hdR['EBV_4102_4341_CORRO2']
			self.hdR['HIERARCH EBV_4102_4341_CORRO3'] = 10**(0.4 *self.hdR['EBV_4102_4341'] *klO3)
			self.hdR['HIERARCH EBV_4102_4341_CORRO3_err']= self.hdR['EBV_4102_4341_err'] * n.log(10) * 0.4 * klO3 * self.hdR['EBV_4102_4341_CORRO3']
			self.hdR['HIERARCH EBV_4102_4341_CORRHb']=10**(0.4 *klHb *self.hdR['EBV_4102_4341'] )
			self.hdR['HIERARCH EBV_4102_4341_CORRHb_err']= self.hdR['EBV_4102_4341_err'] * n.log(10) * 0.4 * klHb * self.hdR['EBV_4102_4341_CORRHb']
		else :
			self.hdR['HIERARCH BD_4102_4341'] = self.dV
			self.hdR['HIERARCH BD_4102_4341_err'] = self.dV
			self.hdR['HIERARCH EBV_4102_4341'] = self.dV
			self.hdR['HIERARCH EBV_4102_4341_err'] = self.dV
			self.hdR['HIERARCH EBV_4102_4341_CORRO2'] = self.dV
			self.hdR['HIERARCH EBV_4102_4341_CORRO2_err'] = self.dV
			self.hdR['HIERARCH EBV_4102_4341_CORRO3'] = self.dV
			self.hdR['HIERARCH EBV_4102_4341_CORRO3_err'] = self.dV
			self.hdR['HIERARCH EBV_4102_4341_CORRHb'] = self.dV
			self.hdR['HIERARCH EBV_4102_4341_CORRHb_err'] = self.dV

		# if BD computation succeeded, we can compute instrinsic quantities
		cor_names = n.array(['EBV_4862_4341_CORR','EBV_4862_4102_CORR', 'EBV_4102_4341_CORR'])
		if len((self.BDarray==1).nonzero()[0]>=1):
			name = cor_names[(self.BDarray==1)][0]
			# intrinsic fluxes O3 4960
			self.hdR['HIERARCH flux_O3_4960_intrinsic'] = self.hdR['O3_4960_flux_nc']/ self.hdR[name+'O3']
			self.hdR['HIERARCH flux_O3_4960_intrinsic_err']= self.hdR['flux_O3_4960_intrinsic'] * ((self.hdR['O3_4960_fluxErr_nc']/ self.hdR['O3_4960_flux_nc'] )**2.+ (self.hdR[name+'O3_err']/ self.hdR[name+'O3'] )**2.)**0.5
			# intrinsic fluxes O3 5007
			self.hdR['HIERARCH flux_O3_5007_intrinsic'] = self.hdR['O3_5007_flux_nc']/ self.hdR[name+'O3']
			self.hdR['HIERARCH flux_O3_5007_intrinsic_err']=self.hdR['flux_O3_5007_intrinsic'] * ((self.hdR['O3_5007_fluxErr_nc']/ self.hdR['O3_5007_flux_nc'] )**2.+ (self.hdR[name+'O3_err']/ self.hdR[name+'O3'] )**2.)**0.5
			# intrinsic fluxes O2
			self.hdR['HIERARCH flux_O2_3728_intrinsic']= self.hdR['O2_3728_flux_nc']/ self.hdR[name+'O2']
			self.hdR['HIERARCH flux_O2_3728_intrinsic_err']=self.hdR['flux_O2_3728_intrinsic'] * ((self.hdR['O2_3728_fluxErr_nc']/ self.hdR['O2_3728_flux_nc'] )**2.+ (self.hdR[name+'O2_err']/ self.hdR[name+'O2'] )**2.)**0.5
			# intrinsic fluxes Hb
			self.hdR['HIERARCH flux_H1_4862_intrinsic'] =self.hdR['H1_4862_flux_nc'] / self.hdR[name+'Hb']
			self.hdR['HIERARCH flux_H1_4862_intrinsic_err'] =self.hdR['flux_H1_4862_intrinsic'] * ((self.hdR['H1_4862_fluxErr_nc']/ self.hdR['H1_4862_flux_nc'])**2.+ (self.hdR[name+'Hb_err'] /self.hdR[name+'Hb'] )**2.)**0.5
			# deduce SFR from [OII] 
			self.hdR['HIERARCH SFR_O2_3728'] = 10**(0.27) * 10**(-41) * self.hdR['flux_O2_3728_intrinsic'] * self.sphereCM.value
			self.hdR['HIERARCH SFR_O2_3728_err'] = self.hdR['SFR_O2_3728'] * self.hdR['flux_O2_3728_intrinsic_err'] / self.hdR['flux_O2_3728_intrinsic']
			# deduce SFR from Hbeta
			self.hdR['HIERARCH SFR_H1_4862'] = 10**(0.58) * 10**(-41) * self.hdR['flux_H1_4862_intrinsic'] * self.sphereCM.value
			self.hdR['HIERARCH SFR_H1_4862_err'] = self.hdR['SFR_H1_4862'] * self.hdR['flux_H1_4862_intrinsic_err'] / self.hdR['flux_H1_4862_intrinsic']
		else :
			self.hdR['HIERARCH flux_O3_4960_intrinsic'] = self.dV
			self.hdR['HIERARCH flux_O3_4960_intrinsic_err'] = self.dV
			self.hdR['HIERARCH flux_O3_5007_intrinsic'] = self.dV
			self.hdR['HIERARCH flux_O3_5007_intrinsic_err'] = self.dV
			self.hdR['HIERARCH flux_O2_3728_intrinsic'] = self.dV
			self.hdR['HIERARCH flux_O2_3728_intrinsic_err'] = self.dV
			self.hdR['HIERARCH flux_H1_4862_intrinsic'] = self.dV
			self.hdR['HIERARCH flux_H1_4862_intrinsic_err'] = self.dV
			self.hdR['HIERARCH SFR_O2_3728'] = self.dV
			self.hdR['HIERARCH SFR_O2_3728_err'] = self.dV
			self.hdR['HIERARCH SFR_H1_4862'] = self.dV
			self.hdR['HIERARCH SFR_H1_4862_err'] = self.dV

		# computes O32
		if self.hdR['O3_4960_flux_nc']>0 and self.hdR['O3_5007_flux_nc']>0 and self.hdR['O2_3728_flux_nc']>0 and self.hdR['O3_4960_fluxErr_nc']>0 and self.hdR['O3_5007_fluxErr_nc']>0 and self.hdR['O2_3728_fluxErr_nc'] >0 :
			self.hdR['HIERARCH O32'] = (self.hdR['O3_4960_flux_nc']+ self.hdR['O3_5007_flux_nc'])/ self.hdR['O2_3728_flux_nc']
			O32ErrFrac =  ( ((self.hdR['O3_4960_fluxErr_nc']+ self.hdR['O3_5007_fluxErr_nc'])/ (self.hdR['O3_4960_flux_nc'] +self.hdR['O3_5007_flux_nc']))**2. + (self.hdR['O2_3728_fluxErr_nc'] /self.hdR['O2_3728_flux_nc'] ) **2.)**0.5  
			self.hdR['HIERARCH O32_err'] = self.hdR['O32'] * O32ErrFrac
		else :
			self.hdR['HIERARCH O32'] = self.dV
			self.hdR['HIERARCH O32_err'] = self.dV

		# compute R23 and 12 log(O/H) with Tremonti 04 estimator
		if self.hdR['O3_4960_flux_nc']>0 and self.hdR['O3_5007_flux_nc']>0 and self.hdR['O2_3728_flux_nc']>0 and self.hdR['H1_4862_flux_nc']>0 and self.hdR['O3_4960_fluxErr_nc']>0 and self.hdR['O3_5007_fluxErr_nc']>0 and self.hdR['O2_3728_fluxErr_nc'] >0 and self.hdR['H1_4862_fluxErr_nc']>0 :
			self.hdR['HIERARCH R23'] = (self.hdR['O3_4960_flux_nc']+self.hdR['O3_5007_flux_nc']+ self.hdR['O2_3728_flux_nc'])/self.hdR['H1_4862_flux_nc']
			R23ErrFrac=( ((self.hdR['O3_4960_fluxErr_nc']+ self.hdR['O3_5007_fluxErr_nc']+ self.hdR['O2_3728_fluxErr_nc']) / (self.hdR['O3_4960_flux_nc']+ self.hdR['O3_5007_flux_nc']+ self.hdR['O2_3728_flux_nc']))**2. + (self.hdR['H1_4862_fluxErr_nc']/ self.hdR['H1_4862_flux_nc'])**2.)**0.5  
			self.hdR['HIERARCH R23_err'] = self.hdR['R23'] * R23ErrFrac
			# 12 log(O/H) with Tremonti 04 estimator
			if self.hdR['R23']>0:
				self.hdR['HIERARCH 12logOH_tremonti04'] = 9.185-0.313*n.log10(self.hdR['R23']) - 0.264 *n.log10(self.hdR['R23'])**2 - 0.321 *n.log10(self.hdR['R23'])**3
				self.hdR['HIERARCH 12logOH_tremonti04_err'] = -0.313* R23ErrFrac / n.log(10) - 0.264 *2 * R23ErrFrac / n.log(10) * n.log10(self.hdR['R23']) - 0.321 * 3* R23ErrFrac / n.log(10) * n.log10(self.hdR['R23'])**2
			else :
				self.hdR['HIERARCH 12logOH_tremonti04'] = self.dV
				self.hdR['HIERARCH 12logOH_tremonti04_err'] = self.dV
				
		else :
			self.hdR['HIERARCH R23'] = self.dV
			self.hdR['HIERARCH R23_err'] = self.dV
			self.hdR['HIERARCH 12logOH_tremonti04'] = self.dV
			self.hdR['HIERARCH 12logOH_tremonti04_err'] = self.dV
		
		"""
		if self.hdR['flux_O3_4960_intrinsic']>0 and self.hdR['flux_O3_5007_intrinsic']>0 and self.hdR['flux_O2_3728_intrinsic']>0 and self.hdR['flux_H1_4862_intrinsic']>0 and self.hdR['flux_O2_3728_intrinsic_err']>0 and self.hdR['flux_O3_5007_intrinsic_err'] >0 :
			# compute R23 intrinsic
			self.hdR['HIERARCH R23_intrinsic'] = (self.hdR['flux_O3_4960_intrinsic']+ self.hdR['flux_O3_5007_intrinsic']+ self.hdR['flux_O2_3728_intrinsic']) /self.hdR['flux_H1_4862_intrinsic']
			R23ErrFrac_intrinsic=( ((self.hdR['flux_O3_4960_intrinsic_err']+ self.hdR['flux_O3_5007_intrinsic_err'] + self.hdR['flux_O2_3728_intrinsic_err']) / (self.hdR['flux_O3_4960_intrinsic'] + self.hdR['flux_O3_5007_intrinsic']+ self.hdR['flux_O2_3728_intrinsic'])) **2. + (self.hdR['flux_H1_4862_intrinsic_err']/ self.hdR['flux_H1_4862_intrinsic']) **2.)**0.5  
			self.hdR['HIERARCH R23_intrinsic_err'] = self.hdR['R23_intrinsic'] * R23ErrFrac_intrinsic
			if self.hdR['R23_intrinsic']>0:
				self.hdR['HIERARCH 12logOH_tremonti04_intrinsic'] = 9.185-0.313* n.log10(self.hdR['R23_intrinsic']) - 0.264 *n.log10(self.hdR['R23_intrinsic'])**2 - 0.321 *n.log10(self.hdR['R23_intrinsic'])**3
				self.hdR['HIERARCH 12logOH_tremonti04_intrinsic_err'] = -0.313* R23ErrFrac_intrinsic / n.log(10) - 0.264 *2 * R23ErrFrac_intrinsic / n.log(10) * n.log10(self.hdR['R23_intrinsic']) - 0.321 * 3* R23ErrFrac_intrinsic / n.log(10) * n.log10(self.hdR['R23_intrinsic'])**2
		else :
			self.hdR['HIERARCH R23_intrinsic'] = self.dV
			self.hdR['HIERARCH R23_intrinsic_err'] = self.dV
			self.hdR['HIERARCH 12logOH_tremonti04_intrinsic'] = self.dV
			self.hdR['HIERARCH 12logOH_tremonti04_intrinsic_err'] = self.dV
		"""
		
		# 12 log OH with O2(3728)/Hbeta
		if  self.hdR['O2_3728_flux_nc']>0 and self.hdR['H1_4862_flux_nc']>0 and self.hdR['O2_3728_fluxErr_nc']>0 and self.hdR['H1_4862_fluxErr_nc']>0:
			OpH = self.hdR['O2_3728_flux_nc']/self.hdR['H1_4862_flux_nc']
			OpHErrFrac = ( (self.hdR['O2_3728_fluxErr_nc']/self.hdR['O2_3728_flux_nc'])**2. + (self.hdR['H1_4862_fluxErr_nc']/self.hdR['H1_4862_flux_nc'])**2.)**0.5  
			OpHErr = OpH * OpHErrFrac
			if OpH>0:
				self.hdR['HIERARCH 12logO2H'] = n.log10(OpH) + 7.637
				self.hdR['HIERARCH 12logO2H_err'] = OpHErrFrac/n.log(10.)
			else :
				self.hdR['HIERARCH 12logO2H'] = self.dV
				self.hdR['HIERARCH 12logO2H_err'] = self.dV
		else :
			self.hdR['HIERARCH 12logO2H'] = self.dV
			self.hdR['HIERARCH 12logO2H_err'] = self.dV

		# 12 log OH with O3(4960+5007)/Hbeta
		if self.hdR['O3_4960_flux_nc']>0 and self.hdR['O3_5007_flux_nc']>0 and self.hdR['H1_4862_flux_nc']>0 and self.hdR['O3_5007_fluxErr_nc'] >0 and self.hdR['O3_4960_fluxErr_nc']>0 and self.hdR['H1_4862_fluxErr_nc']>0:
			O3H = (self.hdR['O3_4960_flux_nc']+ self.hdR['O3_5007_flux_nc'])/ self.hdR['H1_4862_flux_nc']
			O3HErrFrac =  ( ((self.hdR['O3_4960_fluxErr_nc']+self.hdR['O3_5007_fluxErr_nc'])/ (self.hdR['O3_4960_flux_nc']+ self.hdR['O3_5007_flux_nc']))**2. + (self.hdR['H1_4862_fluxErr_nc']/ self.hdR['H1_4862_flux_nc'])**2.)**0.5  
			O3HErr = O3H* O3HErrFrac
			if O3H>0:
				self.hdR['HIERARCH 12logO3H'] = n.log10( O3H )+7.437
				self.hdR['HIERARCH 12logO3H_err'] = O3HErrFrac / n.log(10.)
			else :
				self.hdR['HIERARCH 12logO3H'] = self.dV
				self.hdR['HIERARCH 12logO3H_err'] = self.dV
		else :
			self.hdR['HIERARCH 12logO3H'] = self.dV
			self.hdR['HIERARCH 12logO3H_err'] = self.dV

		# 12 log OH with O3(5007)/Hbeta
		if self.hdR['O3_5007_flux_nc']>0 and self.hdR['H1_4862_flux_nc']>0 and self.hdR['O3_5007_fluxErr_nc'] >0  and self.hdR['H1_4862_fluxErr_nc']>0:
			O35H = (self.hdR['O3_5007_flux_nc'])/self.hdR['H1_4862_flux_nc']
			O35HErrFrac =  ( (self.hdR['O3_5007_fluxErr_nc']/self.hdR['O3_5007_flux_nc'])**2. + (self.hdR['H1_4862_fluxErr_nc']/self.hdR['H1_4862_flux_nc'])**2.)**0.5  
			O35HErr = O35H* O35HErrFrac
			if O35H>0:
				self.hdR['HIERARCH 12logO3_5007_H'] = n.log10( O35H )
				self.hdR['HIERARCH 12logO3_5007_H_err'] = O35HErrFrac / n.log(10.)
			else :
				self.hdR['HIERARCH 12logO3_5007_H'] = self.dV
				self.hdR['HIERARCH 12logO3_5007_H_err'] = self.dV
		else :
			self.hdR['HIERARCH 12logO3_5007_H'] = self.dV
			self.hdR['HIERARCH 12logO3_5007_H_err'] = self.dV

			
	def plot_fit(self):
		"""
		Plots the fit."""

		age ='age = ' +  str(n.round( 10**self.hdu2.header['age_lightW_mean'] ,3))+ '+('+ str(n.round( 10**self.hdu2.header['age_lightW_mean_up']-10**self.hdu2.header['age_lightW_mean'] ,3)) +')-('+str(n.round( 10**self.hdu2.header['age_lightW_mean']-10**self.hdu2.header['age_lightW_mean_low'] ,3))+') Gyr'

		metallicity = 'log(Z/Zsun) = ' + str(n.round( self.hdu2.header['metallicity_lightW_mean'] ,3))+ '+('+ str(n.round( self.hdu2.header['metallicity_lightW_mean_up'] - self.hdu2.header['metallicity_lightW_mean'] ,3)) +')-('+str(n.round( self.hdu2.header['metallicity_lightW_mean'] - self.hdu2.header['metallicity_lightW_mean_low'] ,3))+')'

		mass = 'log(M/Msun) = ' + str(n.round( self.hdu2.header['stellar_mass_mean'] ,3))+ '+('+ str(n.round( self.hdu2.header['stellar_mass_mean_up'] - self.hdu2.header['stellar_mass_mean'] ,3)) +')-('+str(n.round( self.hdu2.header['stellar_mass_mean'] - self.hdu2.header['stellar_mass_mean_low'] ,3))+')'

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
		p.text(0,1,self.baseN+ ", redshift="+str(n.round(self.redshift,3)) )
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
		p.text(0,0,self.baseN)
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


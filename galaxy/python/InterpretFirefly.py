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
cosmo=co.Planck13 #co.FlatLambdaCDM(H0=70,Om0=0.3)
import astropy.units as u
import astropy.io.fits as fits
import numpy as n
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.stats import scoreatpercentile

import astropy.io.fits as fits

"""
from lineListVac import *
allLinesList = n.array([ [Ne3,Ne3_3869,"Ne3_3869","left"], [O3,O3_4363,"O3_4363","right"], [O3,O3_4960,"O3_4960","left"], [O3,O3_5007,"O3_5007","right"], [N2,N2_6549,"N2_6549","left"], [N2,N2_6585,"N2_6585","right"], [H1,H1_3970,"H1_3970","right"], [H1,H1_4102,"H1_4102","right"], [H1,H1_4341,"H1_4341","right"], [H1,H1_4862,"H1_4862","left"], [H1,H1_6564,"H1_6564","left"]]) 
# other lines that are optional
# , [S2,S2_6718,"S2_6718","left"], [S2,S2_6732,"S2_6732","right"], [Ar3,Ar3_7137,"Ar3_7137","left"], [H1,H1_1216,"H1_1216","right"]

doubletList = n.array([[O2_3727,"O2_3727",O2_3729,"O2_3729",O2_mean]])
"""
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
"""
H1=pn.RecAtom('H',1) # Hydrogen Balmer series

bdc0_ref=H1.getEmissivity(1e4, 1e2, lev_i = 4, lev_j = 2) / H1.getEmissivity(1e4, 1e2, lev_i = 3, lev_j = 2)
bdc1_ref=H1.getEmissivity(1e4, 1e2, lev_i = 4, lev_j = 2) / H1.getEmissivity(1e4, 1e2, lev_i = 5, lev_j = 2)
bdc2_ref=H1.getEmissivity(1e4, 1e2, lev_i = 4, lev_j = 2) / H1.getEmissivity(1e4, 1e2, lev_i = 6, lev_j = 2)
bdc3_ref=H1.getEmissivity(1e4, 1e2, lev_i = 4, lev_j = 2) / H1.getEmissivity(1e4, 1e2, lev_i = 7, lev_j = 2)
bdc4_ref=H1.getEmissivity(1e4, 1e2, lev_i = 4, lev_j = 2) / H1.getEmissivity(1e4, 1e2, lev_i = 8, lev_j = 2)
bdc5_ref=H1.getEmissivity(1e4, 1e2, lev_i = 4, lev_j = 2) / H1.getEmissivity(1e4, 1e2, lev_i = 9, lev_j = 2)

bdc23_ref=H1.getEmissivity(1e4, 1e2, lev_i = 5, lev_j = 2)/H1.getEmissivity(1e4, 1e2, lev_i = 6, lev_j = 2)
"""


class InterpretFirefly:
	"""
	This class helps interpret the combination of the fits of the SPM and line fits.

	:param model_file: fits file generated with a LF in a luminosity bin.
	:param cosmo: cosmology class from astropy
	:param dV: default value that hold the place (default : -9999.99) 
	"""
	def __init__(self, model_file, cosmo=cosmo, dV=-9999.99):
		self.model_file = model_file
		self.cosmo = cosmo
		self.dV = dV
		# define global quantities 
		self.baseN = os.path.basename(self.model_file)
		
		# opens the stack
		print " loads the stack :", self.model_file
		hdus = fits.open(self.model_file)
		self.head = hdus[0]
		self.hduStack = hdus[1]
		
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
		headerA =" Redshift age_universe age_lightW_mean age_lightW_err_plus age_lightW_err_minus metallicity_lightW_mean metallicity_lightW_mean_err_plus metallicity_lightW_mean_err_minus age_massW_mean age_massW_err_plus age_massW_err_minus metallicity_massW_mean metallicity_massW_mean_err_plus metallicity_massW_mean_err_minus stellar_mass stellar_mass_err_plus stellar_mass_err_minus spm_EBV nComponentsSSP"
		
		table_entry = [self.redshift, 10**self.hduSPM.header['age_universe'], 10**self.hduSPM.header['age_lightW_mean'], 10**self.hduSPM.header['age_lightW_mean_up']-10**self.hduSPM.header['age_lightW_mean'], 10**self.hduSPM.header['age_lightW_mean']-10**self.hduSPM.header['age_lightW_mean_low'], self.hduSPM.header['metallicity_lightW_mean'], self.hduSPM.header['metallicity_lightW_mean_up'] - self.hduSPM.header['metallicity_lightW_mean'], self.hduSPM.header['metallicity_lightW_mean'] - self.hduSPM.header['metallicity_lightW_mean_low'], 10**self.hduSPM.header['age_massW_mean'], 10**self.hduSPM.header['age_massW_mean_up']-10**self.hduSPM.header['age_massW_mean'], 10**self.hduSPM.header['age_massW_mean']-10**self.hduSPM.header['age_massW_mean_low'], self.hduSPM.header['metallicity_massW_mean'], self.hduSPM.header['metallicity_massW_mean_up'] - self.hduSPM.header['metallicity_massW_mean'], self.hduSPM.header['metallicity_massW_mean'] - self.hduSPM.header['metallicity_massW_mean_low'], self.hduSPM.header['stellar_mass_mean'], self.hduSPM.header['stellar_mass_mean_up'] - self.hduSPM.header['stellar_mass_mean'], self.hduSPM.header['stellar_mass_mean'] - self.hduSPM.header['stellar_mass_mean_low'], self.hduSPM.header['EBV'], self.hduSPM.header['ssp_number']]
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
				table_entry.append([self.dV, self.dV, self.dV, self.dV, self.dV, self.dV])
				headerA += ' stellar_mass_ssp_'+str(iii) + ' age_ssp_'+str(iii) + ' metal_ssp_'+str(iii) + ' SFR_ssp_'+str(iii) + ' weightMass_ssp_'+str(iii) + ' weightLight_ssp_'+str(iii)

		table_entry = n.array( n.hstack((table_entry)) )
		#print table_entry.shape
		return n.hstack((table_entry, self.hduFull.data[0])), headerA
		
	def get_table_entry_line(self):
		headerA =" Redshift age_universe age_lightW_mean age_lightW_err_plus age_lightW_err_minus metallicity_lightW_mean metallicity_lightW_mean_err_plus metallicity_lightW_mean_err_minus age_massW_mean age_massW_err_plus age_massW_err_minus metallicity_massW_mean metallicity_massW_mean_err_plus metallicity_massW_mean_err_minus stellar_mass stellar_mass_err_plus stellar_mass_err_minus spm_EBV nComponentsSSP"
		
		table_entry = [self.redshift, 10**self.hduSPM.header['age_universe'], 10**self.hduSPM.header['age_lightW_mean'], 10**self.hduSPM.header['age_lightW_mean_up']-10**self.hduSPM.header['age_lightW_mean'], 10**self.hduSPM.header['age_lightW_mean']-10**self.hduSPM.header['age_lightW_mean_low'], self.hduSPM.header['metallicity_lightW_mean'], self.hduSPM.header['metallicity_lightW_mean_up'] - self.hduSPM.header['metallicity_lightW_mean'], self.hduSPM.header['metallicity_lightW_mean'] - self.hduSPM.header['metallicity_lightW_mean_low'], 10**self.hduSPM.header['age_massW_mean'], 10**self.hduSPM.header['age_massW_mean_up']-10**self.hduSPM.header['age_massW_mean'], 10**self.hduSPM.header['age_massW_mean']-10**self.hduSPM.header['age_massW_mean_low'], self.hduSPM.header['metallicity_massW_mean'], self.hduSPM.header['metallicity_massW_mean_up'] - self.hduSPM.header['metallicity_massW_mean'], self.hduSPM.header['metallicity_massW_mean'] - self.hduSPM.header['metallicity_massW_mean_low'], self.hduSPM.header['stellar_mass_mean'], self.hduSPM.header['stellar_mass_mean_up'] - self.hduSPM.header['stellar_mass_mean'], self.hduSPM.header['stellar_mass_mean'] - self.hduSPM.header['stellar_mass_mean_low'], self.hduSPM.header['EBV'], self.hduSPM.header['ssp_number']]
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
				table_entry.append([self.dV, self.dV, self.dV, self.dV, self.dV, self.dV])
				headerA += ' stellar_mass_ssp_'+str(iii) + ' age_ssp_'+str(iii) + ' metal_ssp_'+str(iii) + ' SFR_ssp_'+str(iii) + ' weightMass_ssp_'+str(iii) + ' weightLight_ssp_'+str(iii)

		table_entry = n.array( n.hstack((table_entry)) )
		return n.hstack((table_entry, self.hduLine.data[0])), headerA
		
	def print_spm_result(self):
		age ='age = ' +  str(n.round( 10**self.hdu2.header['light_age'] ,3))+ '+('+ str(n.round( 10**self.hdu2.header['age_lightW_mean_up']-10**self.hdu2.header['light_age'] ,3)) +')-('+str(n.round( 10**self.hdu2.header['age_lightW_mean']-10**self.hdu2.header['age_lightW_mean_low'] ,3))+') Gyr'
		metallicity = 'log(Z/Zsun) = ' + str(n.round( self.hdu2.header['metallicity_lightW_mean'] ,3))+ '+('+ str(n.round( self.hdu2.header['metallicity_lightW_mean_up'] - self.hdu2.header['metallicity_lightW_mean'] ,3)) +')-('+str(n.round( self.hdu2.header['metallicity_lightW_mean'] - self.hdu2.header['metallicity_lightW_mean_low'] ,3))+')'
		mass = 'log(M/Msun) = ' + str(n.round( self.hdu2.header['stellar_mass_mean'] ,3))+ '+('+ str(n.round( self.hdu2.header['stellar_mass_mean_up'] - self.hdu2.header['stellar_mass_mean'] ,3)) +')-('+str(n.round( self.hdu2.header['stellar_mass_mean'] - self.hdu2.header['stellar_mass_mean_low'] ,3))+')'
		return age, metallicity, mass
		


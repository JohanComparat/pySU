"""
.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

General purpose:
................

The class StellarPopulationModel is a wrapper dedicated to handling the fit of stellar population models on observed spectra. 
It gathers all inputs : from the model and from the data.

*Imports*::

	import numpy as np
	import astropy.io.fits as pyfits
	import astropy.units as u
	import glob
	import pandas as pd
	import os
	from firefly_instrument import *
	from firefly_dust import *
	from firefly_fitter import *
	from firefly_library import *

"""
import numpy as np
import astropy.io.fits as pyfits
import astropy.units as u
import glob
import pandas as pd
import os
from os.path import join
#from scipy.stats import sigmaclip

#from firefly_dust import *
#import firefly_dust as f_dust
from firefly_dust import hpf, unred, determine_attenuation
from firefly_instrument import downgrade 
from firefly_fitter import fitter
from firefly_library import airtovac, convert_chis_to_probs, light_weights_to_mass, calculate_averages_pdf, normalise_spec, match_data_models


class StellarPopulationModel:
	"""
	:param specObs: specObs observed spectrum object initiated with the  GalaxySpectrumFIREFLY class.
	:param models: choose between 'm11', 'bc03' or 'm09'. 

		* m11 corresponds to all the models compared in `Maraston and Stromback 2011  <http://adsabs.harvard.edu/abs/2011MNRAS.418.2785M>`_. 
		* m09 to `Maraston et al. 2009 <http://adsabs.harvard.edu/abs/2009A%26A...493..425M>`_. 
		* bc03 to the `Bruzual and Charlot 2003 models <http://adsabs.harvard.edu/abs/2003MNRAS.344.1000B>`_.

	:param model_libs: only necessary if using m11. 
	Choose between `MILES <http://adsabs.harvard.edu/abs/2011A%26A...532A..95F>`_, MILES revisednearIRslope, MILES UVextended, `STELIB <http://adsabs.harvard.edu/abs/2003A%26A...402..433L>`_, `ELODIE <http://adsabs.harvard.edu/abs/2007astro.ph..3658P>`_, `MARCS <http://adsabs.harvard.edu/abs/2008A%26A...486..951G>`_. 

		* MILES, MILES revisednearIRslope, MILES UVextended, STELIB, ELODIE are empirical libraries. 
		* MARCS is a theoretical library.

	:param imfs: choose the `initial mass function <https://en.wikipedia.org/wiki/Initial_mass_function>`_:

		* 'ss' for `Salpeter <http://adsabs.harvard.edu/abs/1955ApJ...121..161S>`_or 
		* 'kr' for `Kroupa <http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1112.3340>`_ or 
		* 'cha' for `Chabrier <http://adsabs.harvard.edu/abs/2003PASP..115..763C>`_.

	 Notes
	 -----

	.. note::
		*This is how it proceeds :*
		 #. reads the parameter file by using parameters_obtain(parameters.py)
		 #. It opens the data file, model files, then it matches their resolutions by downgrading the models to instrumental and velocity dispersion resolution
		 #. Determines dust attenuation curve to be applied to the models. Two options : through HPF fitting (3.1.) or through filtered values to determing SP properties (3.2.).
		 #. It fits the models to the data
		 #. Gets mass-weighted SSP contributions using saved M/L ratio.
		 #. Convert chis into probabilities and calculates all average properties and errors (assuming the number of degrees of freedom is the number of wavelength points)
		 #. Optionally produces a plot
		 #. Finally, it writes the output files

	"""
	def __init__(self, specObs, outputFile, cosmo, models = 'm11', model_libs = ['MILES_UVextended'], imfs = ['ss','kr'], hpf_mode = 'on', age_limits = [6,10.1], downgrade_models = True, dust_law = 'calzetti', max_ebv = 1.5, num_dust_vals = 200, dust_smoothing_length = 200, max_iterations = 10, fit_per_iteration_cap = 1000, pdf_sampling = 300, data_wave_medium = 'vacuum', Z_limits = [-0.1,0.1], wave_limits = [0,99999990], suffix = "-fireflyFits.fits",use_downgraded_models = False):
		self.cosmo = cosmo
		self.specObs = specObs
		self.outputFile = outputFile
		#################### STARTS HERE #################### 
		# sets the models
		self.models = models # m11/bc03 / m09
		self.model_libs = model_libs
		self.suffix = suffix
		self.deltal_libs = []
		self.vdisp_round = int(round(self.specObs.vdisp/5.0)*5.0) # rounding vDisp for the models
		self.use_downgraded_models = use_downgraded_models
		if self.models == 'm11':
			for m in self.model_libs:
				if m == 'MILES' or m == 'MILES_revisednearIRslope' or m == 'MILES_UVextended':
					self.deltal_libs.append(2.55)
				elif m == 'STELIB':
					self.deltal_libs.append(3.40)
				elif m == 'ELODIE':
					self.deltal_libs.append(0.55)
				elif m == 'MARCS':
					self.deltal_libs.append(0.1)

		elif self.models=='bc03':
			self.model_libs = ['STELIB_BC03']
			imfs        = ['cha']
			self.deltal_libs = [3.00]

		elif self.models == 'm09':
			self.model_libs = ['M09']
			if downgrade_models:
				self.deltal_libs = [0.4]
			else:
				self.deltal_libs = [3.6]
		# sets the Initial mass function
		self.imfs = imfs
		self.hpf_mode = hpf_mode
		self.age_limits = age_limits

		self.downgrade_models = downgrade_models
		self.dust_law = dust_law
		self.max_ebv = max_ebv
		self.num_dust_vals = num_dust_vals
		self.dust_smoothing_length = dust_smoothing_length
		# Specific fitting options
		self.max_iterations = max_iterations
		self.fit_per_iteration_cap = fit_per_iteration_cap
		# Sampling size when calculating the maximum pdf (100=recommended)
		self.pdf_sampling = pdf_sampling
		# Default is air, unless manga is used
		self.data_wave_medium = data_wave_medium
		self.Z_limits = Z_limits
		self.wave_limits = wave_limits

	def get_model(self, model_used, imf_used, deltal, vdisp, wave_instrument, r_instrument, ebv_mw):

		"""
		Retrieves all relevant model files, in their downgraded format.
		If they aren't downgraded to the correct resolution / velocity dispersion,
		takes the base models in their native form and converts to downgraded files.

		:param model_used: list of models to be used, for example ['m11', 'm09'].
		:param imf_used: list of imf to be used, for example ['ss', 'cha'].
		:param deltal: delta lambda in the models.
		:param vdisp: velocity dispersion observed in the galaxy.
		:param wave_instrument: wavelength array from the observations
		:param r_instrument: resolution array from the observations
		:param  ebv_mw: E(B-V) from the dust maps for the galaxy.
		
		Workflow
		----------		
			A. loads the models m11 or m09: maps parameters to the right files. Then it constructs the model array. Finally converts wavelengths to air or vacuum.
			B. downgrades the model to match data resolution
			C. applies attenuation
			D. stores models in 
				self.model_wavelength, 
				self.model_flux, 
				self.age_model, 
				self.metal_model 
			
			and returns it as well
			
		"""
		# first the m11 case
		if self.models == 'm11':
			first_file  = True 
			model_files = []
			if self.use_downgraded_models :
				if model_used == 'MILES_UVextended' or model_used == 'MILES_revisedIRslope':
					model_path 		= join(os.environ['STELLARPOPMODELS_DIR'],'data','SSP_M11_MILES_downgraded','ssp_M11_' + model_used+ '.' + imf_used)				
				else:
					model_path 		= join(os.environ['STELLARPOPMODELS_DIR'],'data','SSP_M11_'+ model_used + '_downgraded', 'ssp_M11_' +model_used +'.' + imf_used)
			else:
				if model_used == 'MILES_UVextended' or model_used == 'MILES_revisedIRslope':
					model_path 		= join(os.environ['STELLARPOPMODELS_DIR'],'data','SSP_M11_MILES', 'ssp_M11_'+model_used+'.'+imf_used)
				else:
					model_path 		= join(os.environ['STELLARPOPMODELS_DIR'],'data','SSP_M11_'+model_used ,'ssp_M11_' +model_used +'.' + imf_used)

			# Constructs the metallicity array of models :
			all_metal_files = glob.glob(model_path+'*')
			#print all_metal_files
			metal_files 	= []
			metal 	    = []
			for z in range(len(all_metal_files)):
				zchar = all_metal_files[z][len(model_path):]
				if zchar == 'z001':
					znum = -0.3
				elif zchar == 'z002':
					znum = 0.0
				elif zchar == 'z004':
					znum = 0.3
				elif zchar == 'z0001.bhb':
					znum = -1.301
				elif zchar == 'z0001.rhb':
					znum = -1.302
				elif zchar == 'z10m4.bhb':
					znum = -2.301
				elif zchar == 'z10m4.rhb':
					znum = -2.302
				elif zchar == 'z10m4':
					znum = -2.300
				else:
					raise NameError('Unrecognised metallicity! Check model file names.')

				if znum>self.Z_limits[0] and znum<self.Z_limits[1]:
					metal_files.append(all_metal_files[z])
					metal.append(znum)
			# constructs the model array
			model_flux, age_model, metal_model = [],[],[]
			for zi,z in enumerate(metal_files):
				print "Retrieving and downgrading models for "+z
				model_table = pd.read_table(z,converters={'Age':np.float64}, header=None ,usecols=[0,2,3], names=['Age','wavelength_model','flux_model'], delim_whitespace=True)
				age_data = np.unique(model_table['Age'].values.ravel())
				for a in age_data:
					logyrs_a = np.log10(a)+9.0
					#print "age model selection:", self.age_limits[0], logyrs_a, self.age_limits[1]
					if logyrs_a < self.age_limits[0] or logyrs_a > self.age_limits[1]:
						continue

					spectrum = model_table.loc[model_table.Age == a, ['wavelength_model', 'flux_model'] ].values
					wavelength_int,flux = spectrum[:,0],spectrum[:,1]

					# converts to air wavelength
					if self.data_wave_medium == 'vacuum':
						wavelength = airtovac(wavelength_int)
					else:
						wavelength = wavelength_int

					# downgrades the model
					if self.downgrade_models:
						mf = downgrade(wavelength,flux,deltal,self.vdisp_round, wave_instrument, r_instrument)
					else:
						mf = copy.copy(flux)			
			
					# Reddens the models
					if ebv_mw != 0:
						attenuations = unred(wavelength,ebv=0.0-ebv_mw)
						model_flux.append(mf*attenuations)
					else:
						model_flux.append(mf)

					age_model.append(a)
					metal_model.append(metal[zi])
					first_model = False

			print "Retrieved all models!"
			self.model_wavelength, self.model_flux, self.age_model, self.metal_model = wavelength, model_flux, age_model, metal_model
			return wavelength, model_flux, age_model, metal_model

		elif self.models == 'm09':
			first_file  = True 
			model_files = []
			if self.use_downgraded_models:
				model_path = join(os.environ['STELLARPOPMODELS_DIR'],'data', 'UVmodels_Marastonetal08b_downgraded')
			else:
				model_path = join(os.environ['STELLARPOPMODELS_DIR'],'data', 'UVmodels_Marastonetal08b')
			# Gathers the list of models with metallicities and ages of interest:
			all_metal_files = glob.glob(model_path+'*')
			metal_files 	= []
			metal 			= []
			for z in range(len(all_metal_files)):
				zchar = all_metal_files[z].split('.')[1][2:]
				if zchar == 'z001':
					znum = -0.3
				elif zchar == 'z002':
					znum = 0.0
				elif zchar == 'z004':
					znum = 0.3
				elif zchar == 'z0001':
					znum = -1.300
				else:
					raise NameError('Unrecognised metallicity! Check model file names.')

				if znum>self.Z_limits[0] and znum<self.Z_limits[1]:
					metal_files.append(all_metal_files[z])
					metal.append(znum)

			# constructs the model array
			model_flux, age_model, metal_model = [],[],[]
			for zi,z in enumerate(metal_files):
				print "Retrieving and downgrading models for "+z
				model_table = pd.read_table(z,converters={'Age':np.float64}, header=None ,usecols=[0,2,3], names=['Age','wavelength_model','flux_model'], delim_whitespace=True)
				age_data = np.unique(model_table['Age'].values.ravel())
				for a in age_data:
					logyrs_a = np.log10(a)+9.0
					#print "age model selection:", self.age_limits[0], logyrs_a, self.age_limits[1]
					if logyrs_a < self.age_limits[0] or logyrs_a > self.age_limits[1]:
						continue

					spectrum = model_table.loc[model_table.Age == a, ['wavelength_model', 'flux_model'] ].values
					wavelength_int,flux = spectrum[:,0],spectrum[:,1]

					# converts to air wavelength
					if self.data_wave_medium == 'vacuum':
						wavelength = airtovac(wavelength_int)
					else:
						wavelength = wavelength_int

					# downgrades the model
					if self.downgrade_models:
						mf = downgrade(wavelength,flux,deltal,self.vdisp_round, wave_instrument, r_instrument)
					else:
						mf = copy.copy(flux)			
			
					# Reddens the models
					if ebv_mw != 0:
						attenuations = unred(wavelength,ebv=0.0-ebv_mw)
						model_flux.append(mf*attenuations)
					else:
						model_flux.append(mf)

					age_model.append(a)
					metal_model.append(metal[zi])
					first_model = False

			print "Retrieved all models!"
			self.model_wavelength, self.model_flux, self.age_model, self.metal_model = wavelength, model_flux, age_model, metal_model
			return wavelength, model_flux, age_model, metal_model


	def fit_models_to_data(self):
		"""
		Once the data and models are loaded, then execute this function to find the best model. It loops overs the models to be fitted on the data:
		 #. gets the models
		 #. matches the model and data to the same resolution
		 #. normalises the spectra
		"""
		for mi,mm in enumerate(self.model_libs):
			# loop over the models
			for ii in self.imfs:
				# loop over the IMFs
				# A. gets the models
				print "getting the models"
				deltal = self.deltal_libs[mi]
				model_wave_int, model_flux_int, age, metal = self.get_model( mm, ii, deltal, self.specObs.vdisp, self.specObs.restframe_wavelength, self.specObs.r_instrument, self.specObs.ebv_mw)
				# B. matches the model and data to the same resolution
				print "Matching models to data"
				wave, data_flux, error_flux, model_flux_raw = match_data_models( self.specObs.restframe_wavelength, self.specObs.flux, self.specObs.bad_flags, self.specObs.error, model_wave_int, model_flux_int, self.wave_limits[0], self.wave_limits[1], saveDowngradedModel = False)
				# C. normalises the models to the median value of the data
				print "Normalising the models"
				model_flux, mass_factors = normalise_spec(data_flux, model_flux_raw)

			# 3. Corrects from dust attenuation
			if self.hpf_mode=='on':
				# 3.1. Determining attenuation curve through HPF fitting, apply attenuation curve to models and renormalise spectra
				best_ebv, attenuation_curve = determine_attenuation(wave, data_flux, error_flux, model_flux, self, age, metal)
				model_flux_atten = np.zeros(np.shape(model_flux_raw))
				for m in range(len(model_flux_raw)):
					model_flux_atten[m] = attenuation_curve * model_flux_raw[m]

				model_flux, mass_factors = normalise_spec(data_flux, model_flux_atten)
				# 4. Fits the models to the data
				light_weights, chis, branch = fitter(wave, data_flux, error_flux, model_flux, self)
			
			elif self.hpf_mode == 'hpf_only':

				# 3.2. Uses filtered values to determing SP properties only."
				smoothing_length = self.dust_smoothing_length
				hpf_data    = hpf(data_flux)
				hpf_models  = np.zeros(np.shape(model_flux))
				for m in range(len(model_flux)):
					hpf_models[m] = hpf(model_flux[m])

				zero_dat = np.where( (np.isnan(hpf_data)) & (np.isinf(hpf_data)) )
				hpf_data[zero_dat] = 0.0
				for m in range(len(model_flux)):
					hpf_models[m,zero_dat] = 0.0
				hpf_error    = np.zeros(len(error_flux))
				hpf_error[:] = np.median(error_flux)/np.median(data_flux) * np.median(hpf_data)
				hpf_error[zero_dat] = np.max(hpf_error)*999999.9

				best_ebv = -99
				hpf_models,mass_factors = normalise_spec(hpf_data,hpf_models)
				# 4. Fits the models to the data
				light_weights, chis, branch = fitter(wave, hpf_data,hpf_error, hpf_models, self)
			
			# 5. Get mass-weighted SSP contributions using saved M/L ratio.
			unnorm_mass, mass_weights = light_weights_to_mass(light_weights, mass_factors)
			print "Fitting complete"


			print "Calculating average properties and outputting"
			# 6. Convert chis into probabilities and calculates all average properties and errors
			self.dof = len(wave)
			probs = convert_chis_to_probs(chis, self.dof)
			dist_lum	= self.cosmo.luminosity_distance( self.specObs.redshift).to( u.cm ).value
			averages = calculate_averages_pdf(probs, light_weights, mass_weights, unnorm_mass, age, metal, self.pdf_sampling, dist_lum)
			unique_ages 				= np.unique(age)
			marginalised_age_weights 	= np.zeros(np.shape(unique_ages))
			marginalised_age_weights_int = np.sum(mass_weights.T,1)
			for ua in range(len(unique_ages)):
				marginalised_age_weights[ua] = np.sum(marginalised_age_weights_int[np.where(age==unique_ages[ua])])

			best_fit_index = [np.argmin(chis)]
			best_fit = np.dot(light_weights[best_fit_index],model_flux)[0]

			# stores outputs in the object
			self.best_fit_index = best_fit_index 
			self.best_fit = best_fit
			self.model_flux = model_flux
			self.dist_lum = dist_lum
			self.age = np.array(age)
			self.metal = np.array(metal)
			self.mass_weights = mass_weights
			self.light_weights = light_weights
			self.chis = chis
			self.branch = branch
			self.unnorm_mass = unnorm_mass
			self.probs = probs
			self.wave = wave
			self.best_fit = best_fit
			self.averages = averages

			bf_mass = (self.mass_weights[self.best_fit_index]>0)[0]
			bf_light = (self.light_weights[self.best_fit_index]>0)[0]
			mass_per_ssp = self.unnorm_mass[self.best_fit_index[0]][bf_mass]*10.0**(-17) * 4 * np.pi * self.dist_lum**2.0
			age_per_ssp = self.age[bf_mass]*10**9
			metal_per_ssp = self.metal[bf_mass]
			weight_mass_per_ssp = self.mass_weights[self.best_fit_index[0]][bf_mass]
			weight_light_per_ssp = self.light_weights[self.best_fit_index[0]][bf_light]
			order = np.argsort(-weight_light_per_ssp)

			print "M Msun", self.averages['stellar_mass'], np.log10(mass_per_ssp[order])
			print "age Gyr", 10**self.averages['light_age'], 10**self.averages['mass_age'], age_per_ssp[order]/1e9
			print "Z", self.averages['light_metal'], self.averages['mass_metal'], metal_per_ssp[order]
			print "SFR Msun/yr", mass_per_ssp[order]/age_per_ssp[order]
			print "wm", weight_mass_per_ssp[order]
			print "wl", weight_light_per_ssp[order]
			print "z, age Gyr", self.specObs.redshift, self.cosmo.age(self.specObs.redshift).value
			
			# 8. It writes the output file
			waveCol = pyfits.Column(name="wavelength",format="D", unit="Angstrom", array= wave)
			best_fitCol = pyfits.Column(name="firefly_model",format="D", unit="1e-17erg/s/cm2/Angstrom", array= best_fit)
			cols = pyfits.ColDefs([ waveCol, best_fitCol]) 
			tbhdu = pyfits.BinTableHDU.from_columns(cols)
			
			tbhdu.header['HIERARCH age_universe'] = np.log10(self.cosmo.age(self.specObs.redshift).value*10**9)
			tbhdu.header['HIERARCH redshift'] = self.specObs.redshift

			# mean quantities
			tbhdu.header['HIERARCH age_lightW_mean'] = np.log10(10**9 * 10**averages['light_age']) # log(Gyrs)
			tbhdu.header['HIERARCH age_lightW_mean_up'] = np.log10(10**9 * 10**averages['light_age_1_sig_plus']) # log(Gyrs)
			tbhdu.header['HIERARCH age_lightW_mean_low'] = np.log10(10**9 * 10**averages['light_age_1_sig_minus']) # log(Gyrs)
			tbhdu.header['HIERARCH metallicity_lightW_mean'] = averages['light_metal']
			tbhdu.header['HIERARCH metallicity_lightW_mean_up'] = averages['light_metal_1_sig_plus']
			tbhdu.header['HIERARCH metallicity_lightW_mean_low'] = averages['light_metal_1_sig_minus']
			tbhdu.header['HIERARCH age_massW_mean'] = np.log10(10**9 * 10**averages['mass_age']) # log(Gyrs)
			tbhdu.header['HIERARCH age_massW_mean_up'] = np.log10(10**9 * 10**averages['mass_age_1_sig_plus']) # log(Gyrs)
			tbhdu.header['HIERARCH age_massW_mean_low'] = np.log10(10**9 * 10**averages['mass_age_1_sig_minus']) # log(Gyrs)
			tbhdu.header['HIERARCH metallicity_massW_mean'] = averages['mass_metal']
			tbhdu.header['HIERARCH metallicity_massW_mean_up'] = averages['mass_metal_1_sig_plus']
			tbhdu.header['HIERARCH metallicity_massW_mean_low'] = averages['mass_metal_1_sig_minus']
			tbhdu.header['HIERARCH EBV'] = best_ebv
			tbhdu.header['HIERARCH stellar_mass_mean'] = averages['stellar_mass']
			tbhdu.header['HIERARCH stellar_mass_mean_up'] = averages['stellar_mass_1_sig_plus']
			tbhdu.header['HIERARCH stellar_mass_mean_low'] = averages['stellar_mass_1_sig_minus']
			
			tbhdu.header['HIERARCH ssp_number'] =len(order)
			# quantities per SSP
			if len(order)==3:
				tbhdu.header['HIERARCH stellar_mass_ssp_0'] = np.log10(mass_per_ssp[order])[0]
				tbhdu.header['HIERARCH stellar_mass_ssp_1'] = np.log10(mass_per_ssp[order])[1]
				tbhdu.header['HIERARCH stellar_mass_ssp_2'] = np.log10(mass_per_ssp[order])[2]
				tbhdu.header['HIERARCH age_ssp_0'] = np.log10(age_per_ssp[order][0])
				tbhdu.header['HIERARCH age_ssp_1'] = np.log10(age_per_ssp[order][1])
				tbhdu.header['HIERARCH age_ssp_2'] = np.log10(age_per_ssp[order][2])
				tbhdu.header['HIERARCH metal_ssp_0'] = metal_per_ssp[order][0]
				tbhdu.header['HIERARCH metal_ssp_1'] = metal_per_ssp[order][1]
				tbhdu.header['HIERARCH metal_ssp_2'] = metal_per_ssp[order][2]
				tbhdu.header['HIERARCH SFR_ssp_0'] = mass_per_ssp[order][0]/age_per_ssp[order][0]
				tbhdu.header['HIERARCH SFR_ssp_1'] = mass_per_ssp[order][1]/age_per_ssp[order][1]
				tbhdu.header['HIERARCH SFR_ssp_2'] = mass_per_ssp[order][2]/age_per_ssp[order][2]
				tbhdu.header['HIERARCH weightMass_ssp_0'] = weight_mass_per_ssp[order][0]
				tbhdu.header['HIERARCH weightMass_ssp_1'] = weight_mass_per_ssp[order][1]
				tbhdu.header['HIERARCH weightMass_ssp_2'] = weight_mass_per_ssp[order][2]
				tbhdu.header['HIERARCH weightLight_ssp_0'] = weight_light_per_ssp[order][0]
				tbhdu.header['HIERARCH weightLight_ssp_1'] = weight_light_per_ssp[order][1]
				tbhdu.header['HIERARCH weightLight_ssp_2'] = weight_light_per_ssp[order][2]

			if len(order)==2:
				tbhdu.header['HIERARCH stellar_mass_ssp_0'] = np.log10(mass_per_ssp[order])[0]
				tbhdu.header['HIERARCH stellar_mass_ssp_1'] = np.log10(mass_per_ssp[order])[1]
				tbhdu.header['HIERARCH stellar_mass_ssp_2'] = 0.
				tbhdu.header['HIERARCH age_ssp_0'] = np.log10(age_per_ssp[order][0])
				tbhdu.header['HIERARCH age_ssp_1'] = np.log10(age_per_ssp[order][1])
				tbhdu.header['HIERARCH age_ssp_2'] = 0.
				tbhdu.header['HIERARCH metal_ssp_0'] = metal_per_ssp[order][0]
				tbhdu.header['HIERARCH metal_ssp_1'] = metal_per_ssp[order][1]
				tbhdu.header['HIERARCH metal_ssp_2'] = 0.
				tbhdu.header['HIERARCH SFR_ssp_0'] = mass_per_ssp[order][0]/age_per_ssp[order][0]
				tbhdu.header['HIERARCH SFR_ssp_1'] = mass_per_ssp[order][1]/age_per_ssp[order][1]
				tbhdu.header['HIERARCH SFR_ssp_2'] = 0.
				tbhdu.header['HIERARCH weightMass_ssp_0'] = weight_mass_per_ssp[order][0]
				tbhdu.header['HIERARCH weightMass_ssp_1'] = weight_mass_per_ssp[order][1]
				tbhdu.header['HIERARCH weightMass_ssp_2'] = 0.
				tbhdu.header['HIERARCH weightLight_ssp_0'] = weight_light_per_ssp[order][0]
				tbhdu.header['HIERARCH weightLight_ssp_1'] = weight_light_per_ssp[order][1]
				tbhdu.header['HIERARCH weightLight_ssp_2'] = 0.
				
			if len(order)==1:
				tbhdu.header['HIERARCH stellar_mass_ssp_0'] = np.log10(mass_per_ssp[order])[0]
				tbhdu.header['HIERARCH stellar_mass_ssp_1'] = 0.
				tbhdu.header['HIERARCH stellar_mass_ssp_2'] = 0.
				tbhdu.header['HIERARCH age_ssp_0'] = np.log10(age_per_ssp[order][0])
				tbhdu.header['HIERARCH age_ssp_1'] = 0.
				tbhdu.header['HIERARCH age_ssp_2'] = 0.
				tbhdu.header['HIERARCH metal_ssp_0'] = metal_per_ssp[order][0]
				tbhdu.header['HIERARCH metal_ssp_1'] = 0.
				tbhdu.header['HIERARCH metal_ssp_2'] = 0.
				tbhdu.header['HIERARCH SFR_ssp_0'] = mass_per_ssp[order][0]/age_per_ssp[order][0]
				tbhdu.header['HIERARCH SFR_ssp_1'] = 0.
				tbhdu.header['HIERARCH SFR_ssp_2'] = 0.
				tbhdu.header['HIERARCH weightMass_ssp_0'] = weight_mass_per_ssp[order][0]
				tbhdu.header['HIERARCH weightMass_ssp_1'] = 0.
				tbhdu.header['HIERARCH weightMass_ssp_2'] = 0.
				tbhdu.header['HIERARCH weightLight_ssp_0'] = weight_light_per_ssp[order][0]
				tbhdu.header['HIERARCH weightLight_ssp_1'] = 0.
				tbhdu.header['HIERARCH weightLight_ssp_2'] = 0.

			prihdr = pyfits.Header()
			prihdr['file'] = self.specObs.path_to_spectrum
			prihdr['model'] = self.models
			prihdr['ageMin'] = self.age_limits[0]
			prihdr['ageMax'] = self.age_limits[1]
			prihdr['Zmin'] = self.Z_limits[0]
			prihdr['Zmax'] = self.Z_limits[1]
			prihdu = pyfits.PrimaryHDU(header=prihdr)

			thdulist = pyfits.HDUList([prihdu, tbhdu])
			if os.path.isfile(self.outputFile + self.suffix ):
				os.remove(self.outputFile + self.suffix )
			#print self.outputFile + self.suffix , thdulist, thdulist[1].data, thdulist[0].header
			thdulist.writeto(self.outputFile + self.suffix )

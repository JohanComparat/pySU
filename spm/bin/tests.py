#! /usr/bin/env python
from os.path import join
import os
import glob
import numpy as np
import astropy.cosmology as co
cosmo = co.Planck15
#cosmo=co.FlatLambdaCDM(H0=70,Om0=0.3)

# for one galaxy spectrum
import GalaxySpectrumFIREFLY as gs
import StellarPopulationModel as spm

spec_deep2_file = "D:\data\spectraStacks\content\data\O2_3728\O2_3728-DEEP2R24.2-z0.927_stack_N_100_R_8686_L_1.84133041845e+41.fits"
spec_vvds_file = "D:\data\spectraStacks\content\data\O2_3728\O2_3728-VVDSWIDEI22.5-z0.638_stack_N_100_R_8686_L_2.52410366696e+41.fits"


spec_deep2=gs.GalaxySpectrumFIREFLY(spec_deep2_file, milky_way_reddening=False, stack_resolution = 8686. )
spec_deep2.openObservedStack()
file =os.path.basename(spec_deep2_file)
outFile = join(os.environ['SPECTRASTACKS_DIR'], "fits", file[:7], file)
model_deep2 = spm.StellarPopulationModel(spec_deep2, outFile, cosmo, models = 'm11', model_libs = ['MILES'], imfs = ['ss'], age_limits = [6,10], downgrade_models = True, data_wave_medium = 'air', Z_limits = [-3.,1.],suffix="-SPM-MILES.fits", use_downgraded_models = False)
model_deep2.fit_models_to_data()


spec_vvds=gs.GalaxySpectrumFIREFLY(spec_vvds_file, milky_way_reddening=False, stack_resolution = 600. )
spec_vvds.openObservedStack()
file =os.path.basename(spec_vvds_file)
outFile = join(os.environ['SPECTRASTACKS_DIR'], "fits", file[:7], file)
model_vvds = spm.StellarPopulationModel(spec_vvds, outFile, cosmo, models = 'm11', model_libs = ['MILES'], imfs = ['ss'], age_limits = [6,10], downgrade_models = True, data_wave_medium = 'air', Z_limits = [-3.,1.],suffix="-SPM-MILES.fits", use_downgraded_models = False)
model_vvds.fit_models_to_data()

spec_cmass_file = "D:\data\spectraStacks\sdss\spec-3586-55181-0007.fits"
spec_cmass=gs.GalaxySpectrumFIREFLY(spec_cmass_file, milky_way_reddening=True)
spec_cmass.openObservedSDSSSpectrum()
outFile = "D:\data\spectraStacks\sdss\spec-3586-55181-0007.model"
model_sdss = spm.StellarPopulationModel(spec_cmass, outFile, cosmo, models = 'm11', model_libs = ['MILES'], imfs = ['ss'], age_limits = [6,10], downgrade_models = True, data_wave_medium = 'air', Z_limits = [-3.,1.],suffix="-SPM-MILES.fits", use_downgraded_models = False)
model_sdss.fit_models_to_data()


spec_elg_file = "D:\data\spectraStacks\sdss\spec-4220-55447-0209.fits"
spec_elg=gs.GalaxySpectrumFIREFLY(spec_elg_file, milky_way_reddening=True)
spec_elg.openObservedSDSSSpectrum()
outFile = "D:\data\spectraStacks\sdss\spec-4220-55447-0209.model"
model_elg = spm.StellarPopulationModel(spec_elg, outFile, cosmo, models = 'm11', model_libs = ['MILES'], imfs = ['ss'], age_limits = [6,10], downgrade_models = True, data_wave_medium = 'air', Z_limits = [-3.,1.],suffix="-SPM-MILES.fits", use_downgraded_models = False)
model_elg.fit_models_to_data()

Welcome to the stellar population module 
======================================================================

Modules and exports if you run on sciama, ICG portsmouth.
-------------------------------------------------------------------

module load apps/anaconda/2.4.1
module load apps/python/2.7.8/gcc-4.4.7

# pathes to the pySU version you have locally
export PYTHONPATH=$PYTHONPATH:/users/comparat/pySU/galaxy/python/
export PYTHONPATH=$PYTHONPATH:/users/comparat/pySU/spm/python/
export PYSU_DIR='/users/comparat/pySU'
export PYSU_GAL_DIR='/users/comparat/pySU/galaxy'
export PYSU_SPM_DIR='/users/comparat/pySU/spm'

# SPM DATA EXPORTS: where the models are
export STELLARPOPMODELS_DIR='/users/comparat/stellarpopmodels/trunk'

# DATA STRUCTURE: whre the data is
export DATA_DIR='/mnt/lustre/sdss-dr12'

export VVDS_DIR='/mnt/lustre/sdss-dr12/vvds'
export VIPERS_DIR='/mnt/lustre/sdss-dr12/vipers'
export DEEP2_DIR='/mnt/lustre/sdss-dr12/deep2'

export EBOSSDR14_DIR='/mnt/lustre/sdss-dr12/eBOSS-DR14'
export BOSSDR12_DIR='/mnt/lustre/sdss-dr12/BOSS-DR12'
export SDSSDR12_DIR='/mnt/lustre/sdss-dr12/SDSS-DR12'


Organisation of the DATA 
============================

in $DATA_DIR, there is a directory per survey, say $SURVEY. 
Currently, we have VVDS, VIPERS, DEEP2, SDSS DR12, BOSS DR12 and eBOSS DR14.

in each $SURVEY, you have
 - spectra               : where the spectroscopic data is
 - catalogs              : summary catalogs 
 - stellarpop-$RUN  : outputs of the stellarpopulation models
 
In each stellarpop-$RUN dir, you have
 - stellarpop            : model spectrum of the stellar continuum, information at the individual spectrum level
 - catalogs               : catalogs of properties per larger chunks of data, say about 1000 spectra (all columns are kept)
 - flyAll_catalogs      : few summary catalog (less columns)
 - emission_model    : model spectrum of the gas-phase spectrum, once the stellar component is removed.
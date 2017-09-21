#!/bin/bash

# on ds52, cd ~/software/linux/pySU/spm/bin_SMF

# give overall numbers, table 1 and 2
# how many spectra were observed
# Catalog &
# spectrograph & 
# data release & 
# All spectra & 
# Galaxy spectra & 
# with SPM 
python create_table1.py
# writes here os.environ['OBS_REPO'], 'spm', 'results', "table_1.tex" and "table_2.tex"

# Figure 1 
python object_types_mass.py

# figure 2 is an example spectrum created by firefly code

# Creates figure 3 about sigma_m and SNR per pixel
python pdf_SM_error.py

# estimate the mass density probed by SDSS and BOSS and DEEP2

sh run_smf_plots.sh
python smf_plot.py 0.78 0.83 41.8
python smf_sdss_eboss.py 0.0 1.7



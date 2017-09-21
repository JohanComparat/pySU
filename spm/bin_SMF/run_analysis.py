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

# figure 1 is an example spectrum created by firefly code

# Creates figure 2 about sigma_m and SNR per pixel
python pdf_SM_error.py

# Figure 3 
plot M vs Z for all accurate stellar masses
 - SDSS
 - DEEP2
 - COSMOS: or not

z_bins = n.arange(0, 1.7, 0.05)
m_bins = n.arange(7, 13, 0.1)
weight = area 
area_boss = 10000.
area_sdss = 7900.
area_deep2 = xxx
n.histogram2d(redshift, mass, weight, z_bins, m_bins)

density_bins_per_deg2 = n.array([0.1, 1, 10, 100, 1000])

 
N(z) for all accurate stellar masses

sh run_smf_plots.sh
python smf_plot.py 0.78 0.83 41.8
python smf_sdss_eboss.py 0.0 1.7



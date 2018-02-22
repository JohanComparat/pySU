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
python create_table_completeness.py sdss
python create_table_completeness.py boss
python create_table_snr.py sdss
python create_table_snr.py boss
python create_table_delta_mag.py sdss
python create_table_delta_mag.py boss

rm ~/wwwDir/firefly_data/dr14/v1_1_0/tables/*.tex
cp /data36s/comparat/spm/results/*.tex ~/wwwDir/firefly_data/dr14/v1_1_0/tables/

python create_table_completeness.py deep2

# writes here os.environ['OBS_REPO'], 'spm', 'results', "table_1.tex" and "table_2.tex" and *.tex for the appendix tables.

# to be executed only once 
# python measure_SNMEDIAN_DEEP2.py

# Figure 1 
# first row
python3.4 object_types_mass.py
# second row
python3.4 object_types_SNMEDIANALL.py

rm /home/comparat/wwwDir/firefly_data/dr14/v1_1_0/plots/*.png
cp /data36s/comparat/spm/results/mass-redshift-presentation/*.png /home/comparat/wwwDir/firefly_data/dr14/v1_1_0/plots/

# figure 2 is an example spectrum created by firefly code

# Creates figure 3 about sigma_m and SNR per pixel
python pdf_SM_error.py

# Figure 4, stellar mass functions probed by SDSS and BOSS and DEEP2 in different redshift bins
python mass_density.py 0.2 0.5
python mass_density.py 0.5 0.8
python mass_density.py 0.8 1.1

# Figure 5 stellar mass functions probed by [OII] emitters in DEEP2.
python smf_plot.py 0.78 0.83 41.8



#!/bin/bash

# on ds52, 

cd /home/comparat/software/linux/pySU/spm/bin_SMF

python make_fly_plate_list.py v5_10_0
python summary_table_concatenate_spFlyAll.py v5_10_0
python summary_table_merge_specObjAll_SNRall.py v5_10_0
python summary_table_merge_firefly_MAG.py v5_10_0

python make_fly_plate_list.py 26
python summary_table_concatenate_spFlyAll.py 26
python summary_table_merge_specObjAll_SNRall.py 26
python summary_table_merge_firefly_MAG.py 26

#CATALOG are written here 
#/data37s/SDSS/26/catalogs
#/data37s/SDSS/v5_10_0/catalogs

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

cp /data36s/comparat/spm/results/*.tex /data42s/comparat/firefly/v1_1_0/tables/

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
cp /data36s/comparat/spm/results/mass-snr/*.png /home/comparat/wwwDir/firefly_data/dr14/v1_1_0/plots/

# figure 2 is an example spectrum created by firefly code

# Creates figure 3 about sigma_m and SNR per pixel
python pdf_SM_error.py

# Figure 4, stellar mass functions probed by SDSS and BOSS and DEEP2 in different redshift bins
python mass_density.py 0.2 0.5
python mass_density.py 0.5 0.8
python mass_density.py 0.8 1.1

# Figure 5 stellar mass functions probed by [OII] emitters in DEEP2.
python smf_plot.py 0.78 0.83 41.8



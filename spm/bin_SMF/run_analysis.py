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

sh run_smf_plots.sh
python smf_plot.py 0.78 0.83 41.8
python smf_sdss_eboss.py 0.0 1.7



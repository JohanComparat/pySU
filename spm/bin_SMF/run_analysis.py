#!/bin/bash

# give overall numbers, table 1 and 2
# how many spectra were observed

# data to write in the Table:
# Catalog &
# spectrograph & 
# data release & 
# All spectra & 
# Galaxy spectra & 
# with SPM 
python create_table1.py

sh run_smf_plots.sh
python smf_plot.py 0.78 0.83 41.8
python smf_sdss_eboss.py 0.0 1.7



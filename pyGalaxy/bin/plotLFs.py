"""
This script produces quality plots to check that the LFs are fine.
"""

import sys
from lib_plot import *
from lineListAir import *
SNlim = 5


plotDir="/home/comparat/database/VVDS/products/emissionLineLuminosityFunctions/plotLFs-raw/"
dir="/home/comparat/database/VVDS/products/emissionLineLuminosityFunctions/??_????/"
import glob
lf_measurement_files=glob.glob(dir+"*.txt")

for jj in range(len(lf_measurement_files)):
	lf_fits_file=lf_measurement_files[jj][:-4]+".fits"
	print lf_fits_file
	lf_measurement_file=lf_measurement_files[jj]
	plot_EW_LF_measurement(lf_fits_file,lf_measurement_file,plotDir)

plotDir="/home/comparat/database/DEEP2/products/emissionLineLuminosityFunctions/plotLFs-raw/"
dir="/home/comparat/database/DEEP2/products/emissionLineLuminosityFunctions/??_????/"
import glob
lf_measurement_files=glob.glob(dir+"*.txt")

for jj in range(len(lf_measurement_files)):
	lf_fits_file=lf_measurement_files[jj][:-4]+".fits"
	print lf_fits_file
	lf_measurement_file=lf_measurement_files[jj]
	plot_EW_LF_measurement(lf_fits_file,lf_measurement_file,plotDir)

plotDir="/home/comparat/database/VIPERS/products/emissionLineLuminosityFunctions/plotLFs-raw/"
dir="/home/comparat/database/VIPERS/products/emissionLineLuminosityFunctions/??_????/"
import glob
lf_measurement_files=glob.glob(dir+"*.txt")

for jj in range(len(lf_measurement_files)):
	lf_fits_file=lf_measurement_files[jj][:-4]+".fits"
	print lf_fits_file
	lf_measurement_file=lf_measurement_files[jj]
	plot_EW_LF_measurement(lf_fits_file,lf_measurement_file,plotDir)

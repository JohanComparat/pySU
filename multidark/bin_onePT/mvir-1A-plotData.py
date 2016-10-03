from os.path import join
import numpy as n
import astropy.io.fits as fits
import os
import sys

import lib_functions_1pt as lib

import astropy.cosmology as co
cosmo = co.Planck13
import astropy.units as uu

import matplotlib
matplotlib.use('pdf')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p

#Quantity studied
qty = "mvir"
# working directory
dir = join(os.environ['MULTIDARK_LIGHTCONE_DIR'], qty)
# loads summary file
data = fits.open( join(dir, "MD_"+qty+"_summary.fits"))[1].data

NminCount = 1000
logNpmin = 3

zmin = -0.01
zmax = 2.5


# x coordinates definition
logsig = -n.log10(data['sigmaM'])#
lognu = n.log10(data['nu2']**0.5)
#log_mvir = data["log_"+qty]
log_mvir = data["log_"+qty] - n.log10(cosmo.h)
mvir = 10**data["log_"+qty] / cosmo.h

#=======================
#=======================
cos = 'cen'
#=======================
#=======================
# y coordinates
ff = mvir *  data["dNdlnM_"+cos] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
ff_c = mvir *  data["dNdlnM_"+cos+"_c"] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
log_MF = n.log10( ff )
log_MF_c = n.log10(  ff_c )

# redshift selection
zSel = lib.zSelection( data, zmin, zmax )
# mass selection
mSel = lib.mSelection(data, qty, logNpmin)
# minimum number counts selection
nSelCen = lib.nSelection(data, NminCount, cos )
# altogether
ok = (zSel) & (mSel) & (nSelCen)
# selection per box :
MD04=(data["boxName"]=='MD_0.4Gpc')
MD10=(data["boxName"]=='MD_1Gpc_new_rockS')
MD25=(data["boxName"]=='MD_2.5Gpc')
MD40=(data["boxName"]=='MD_4Gpc')
MD25NW=(data["boxName"]=='MD_2.5GpcNW')
MD40NW=(data["boxName"]=='MD_4GpcNW')

# NOW PLOTTING ALL THE DATA
lib.plot_mvir_function_data(log_mvir[ok], logsig[ok], lognu[ok], log_MF[ok], log_MF_c[ok], data['redshift'][ok], zmin, zmax, cos = cos)

lib.plot_mvir_function_data_perBox(log_mvir, log_MF, MD04, MD10, MD25, MD25NW, MD40, MD40NW, cos=cos)

# ERROR PLOT: JK vs. POISSON
x = data["std90_pc_"+cos] 
y = data["dN_counts_"+cos]**(-0.5)
lib.plot_jackknife_poisson_error(x, y, MD04, MD10, MD25, MD25NW, MD40, MD40NW, cos = cos, dir=join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'mvir'))

#=======================
#=======================
cos = 'sat'
#=======================
#=======================
# y coordinates
ff = mvir *  data["dNdlnM_"+cos] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
ff_c = mvir *  data["dNdlnM_"+cos+"_c"] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
log_MF = n.log10( ff )
log_MF_c = n.log10(  ff_c )

# minimum number counts selection
nSelSat = lib.nSelection(data, NminCount, cos )
# altogether
ok = (zSel) & (mSel) & (nSelSat)

# NOW PLOTTING ALL THE DATA
lib.plot_mvir_function_data(log_mvir[ok], logsig[ok], lognu[ok], log_MF[ok], log_MF_c[ok], data['redshift'][ok], zmin, zmax, cos = cos)

lib.plot_mvir_function_data_perBox(log_mvir, log_MF, MD04, MD10, MD25, MD25NW, MD40, MD40NW, cos=cos)

# ERROR PLOT: JK vs. POISSON
x = data["std90_pc_"+cos] 
y = data["dN_counts_"+cos]**(-0.5)
lib.plot_jackknife_poisson_error(x, y, MD04, MD10, MD25, MD25NW, MD40, MD40NW, cos = cos, dir=join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'mvir'))


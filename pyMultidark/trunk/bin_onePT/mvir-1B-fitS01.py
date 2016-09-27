from os.path import join
import numpy as n
import astropy.io.fits as fits
import os
import lib_functions_1pt as lib

import astropy.cosmology as co
cosmo = co.Planck13

import matplotlib
#matplotlib.use('pdf')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p

from scipy.interpolate import interp1d
from scipy.misc import derivative
import astropy.units as uu


sigma = n.arange(0.05,10,0.05)
delta_c = 1.686
f_ST01 = lambda sigma, A, a, p: A * ((2. * a * (delta_c/sigma)**2.) / (  n.pi))**(0.5) * ( 1 + (a*(delta_c/sigma)**2.) **(-p) ) * n.e**( - a * (delta_c/sigma)**2. / 2.)
log_f_ST01 = lambda logSigma, p : n.log10( f_ST01(10.**logSigma, p[0], p[1], p[2]) )

p0 = [0.287, 0.903, 0.322]

#Quantity studied
qty = "mvir"
# working directory
dir = join(os.environ['MULTIDARK_LIGHTCONE_DIR'], qty)
# loads summary file
data = fits.open( join(dir, "MD_"+qty+"_summary.fits"))[1].data

NminCount = 100 # 1000
Npmin =  300 # 1
nolim = [0,1e17]
limits_04 =  n.log10([Npmin*9.63 * 10**7, 5e12])
limits_10 =  n.log10([Npmin*1.51 * 10**9., 5e13])
limits_25 =  n.log10([Npmin*2.359 * 10**10., 5e14])
limits_40 =  n.log10([Npmin* 9.6 * 10**10. , 5e15])
#MPART = n.array([9.63 * 10**7, 1.51 * 10**9, 2.359 * 10**10, 9.6 * 10**10])

zmin = -0.01
zmax = 0.001

#=======================
#=======================
cos = 'cen'
#=======================
#=======================

# redshift selection
zSel = lib.zSelection( data, zmin, zmax )
# mass selection
mSel = lib.mSelection(data, qty, limits_04, limits_10, limits_25,limits_40) 
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

# x coordinates definition
logsig = n.log10(data['sigmaM'])#
log_mvir = (data["log_"+qty+"_min"]+data["log_"+qty+"_max"])/2.
mvir = 10**log_mvir

# y coordinates
rhom = cosmo.critical_density(data["redshift"]).to(uu.solMass/(uu.Mpc)**3.)/(cosmo.H(data["redshift"])/(100*uu.km/(uu.Mpc*uu.s)))**1.
log_MF = n.log10( mvir * data["dNdVdlnM_"+cos]/ rhom.value )
log_MF_c = n.log10(  data["dNdVdlnM_"+cos+"_c"])
log_f =  n.log10(mvir * data["dNdVdlnM_"+cos]/ rhom.value  / abs(data["dlnsigmaM1_o_dlnM"]))
log_f_c =  n.log10(mvir * data["dNdVdlnM_"+cos+"_c"]/ rhom.value  / abs(data["dlnsigmaM1_o_dlnM"]))

# error on y position
#=================
error = data["dN_counts_"+cos]**(-0.5)


pOpt, pCov = lib.fit_mvir_function_z0(data[ok], x_data = logsig[ok], y_data = log_MF[ok], y_err = error[ok], p0 = p0, cos = cos, mode = "curve_fit")


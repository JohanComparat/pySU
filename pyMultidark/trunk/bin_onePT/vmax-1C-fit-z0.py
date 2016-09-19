from os.path import join
import numpy as n
import astropy.io.fits as fits
import os

import astropy.cosmology as co
cosmo = co.Planck13

import matplotlib
#matplotlib.use('pdf')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p

from scipy.optimize import minimize
from scipy.optimize import curve_fit

from scipy.interpolate import interp1d
from scipy.misc import derivative
# working directory
#=================
dir = join(os.environ['MULTIDARK_LIGHTCONE_DIR'], qty)

#Quantity studied
#=================
qty = "vmax"

# fitting function parameters
#=================
NminCount = 1000
limits_04 = [100, 400]
limits_10 = [250, 1000]
limits_25 = [600, 1300]
limits_40 = [1200, 1600]

p0 = n.array([-3, 3., 0.3, 1.])
errorFactor = 3.
systError = 0.01

zmin = -0.01
zmax = 0.001

#=================
# DATA
#=================
# loads summary file
data = fits.open( join(dir, "MD_"+qty+"_summary.fits"))[1].data
# redshift selection
zSel = lib.zSelection( data, zmin, zmax )
# mass selection
mSel = lib.mSelection(data, limits_04, limits_10, limits_25,limits_40) 
# minimum number counts selection
nSel = lib.nSelection(data, NminCount )
# altogether
ok = (zSel) & (mSel) & (nSel)


#=================
#=================
cos = 'cen'
#=================
#=================

# x coordinates definition
#=================
log_vmax = (n.log10(data["log_"+qty+"_min"])+n.log10(data["log_"+qty+"_max"]))/2.
vmax = 10**log_vmax
#print len(vmax), n.min(vmax), n.max(vmax)

# y coordinates
#=================
norm = (100)**3. /(cosmo.H(data["redshift"]).value)**6.
log_VF = n.log10( norm * vmax**3. * data["dNdVdlnM_"+cos])
log_VF_c = n.log10( norm * vmax**3. * data["dNdVdlnM_"+cos+"_c"])
#print n.min(log_VF), n.max(log_VF)

# error on y position
#=================
error = data["dN_counts_"+cos]**(-0.5)
#error_04 = (data['std90_pc_cen'][MD04]*errorFactor+systError)/ n.log(10.) 

pOpt, pCov = lib.fit_vmax_function_z0(data[ok], x_data = log_vmax[ok], y_data = log_VF[ok], y_err = error[ok], p0 = p0, cos = cos, mode = "curve_fit")


#=================
#=================
cos = 'sat'
#=================
#=================

# x coordinates definition
#=================
log_vmax = (n.log10(data["log_"+qty+"_min"])+n.log10(data["log_"+qty+"_max"]))/2.
vmax = 10**log_vmax
#print len(vmax), n.min(vmax), n.max(vmax)


# y coordinates
#=================
norm = (100)**3. /(cosmo.H(data["redshift"]).value)**6.
log_VF = n.log10( norm * vmax**3. * data["dNdVdlnM_"+cos])
log_VF_c = n.log10( norm * vmax**3. * data["dNdVdlnM_"+cos+"_c"])
#print n.min(log_VF), n.max(log_VF)

# error on y position
#=================
error = data["dN_counts_"+cos]**(-0.5)
#error_04 = (data['std90_pc_cen'][MD04]*errorFactor+systError)/ n.log(10.) 

pOpt, pCov = lib.fit_vmax_function_z0(data[ok], x_data = log_vmax[ok], y_data = log_VF[ok], y_err = error[ok], p0 = p0, cos = cos, mode = "curve_fit")








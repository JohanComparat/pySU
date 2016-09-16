import glob
import sys
import cPickle
from os.path import join
import numpy as n
import astropy.io.fits as fits
import os

import astropy.cosmology as co
cosmo = co.Planck13
import astropy.units as uu

import matplotlib
#matplotlib.use('pdf')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p

from scipy.optimize import minimize
from scipy.optimize import curve_fit

from scipy.interpolate import interp1d
from scipy.misc import derivative

dir='..'
dir_04 = join(dir,"MD_0.4Gpc")
dir_10 = join(dir,"MD_1Gpc")
dir_25 = join(dir,"MD_2.5Gpc")
dir_40 = join(dir,"MD_4Gpc")
dir_25N = join(dir,"MD_2.5GpcNW")
dir_40N = join(dir,"MD_4GpcNW")

data = fits.open( join("..", "vmax", "MD_vmax_summary.fits") )[1].data

NminCount = 1000
limits_04 = [100, 400]
limits_10 = [250, 1000]
limits_25 = [600, 1300]
limits_40 = [1200, 1600]
names = n.array(["SMD", "MDPL", "BigMD", "HMD", "BigMDNW", "HMDNW"])

zmin = -0.01
zmax = 0.001
qty = 'vmax'
cos = "cen"

# redshift selection
zSel = (data["redshift"]>zmin)&(data["redshift"]<zmax)
# mass selection
if  cos == "cen":
	mSel = ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>limits_04[0]) &(data["log_"+qty+"_max"]<limits_04[1])) | ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>limits_10[0]) &(data["log_"+qty+"_max"]<limits_10[1])) |  ((data["boxLength"]==2500.)&(data["log_"+qty+"_min"]>limits_25[0]) &(data["log_"+qty+"_max"]<limits_25[1])) |  ((data["boxLength"]==4000.)&(data["log_"+qty+"_min"]>limits_40[0])&(data["log_"+qty+"_max"]<limits_40[1])) 
if  cos == "sat":
	mSel = ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>limits_04[0]) &(data["log_"+qty+"_max"]<limits_04[1])) | ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>limits_10[0]) &(data["log_"+qty+"_max"]<limits_10[1])) 

# minimum number counts selection
nSel = (data['dN_counts_'+cos]>NminCount)
# altogether
ok = (zSel) & (mSel) & (nSel)

# x coordinates
log_vmax = (n.log10(data["log_"+qty+"_min"][ok])+n.log10(data["log_"+qty+"_max"][ok]))/2.
vmax = 10**log_vmax
print len(vmax), n.min(vmax), n.max(vmax)
# y coordinates
log_VF = n.log10( vmax**4. * data["dNdVdlnM_"+cos][ok])
log_VF_c = n.log10( vmax**4. * data["dNdVdlnM_"+cos+"_c"][ok])
print n.min(log_VF), n.max(log_VF)

MD04=(data["boxLength"]==400.)
MD10=(data["boxLength"]==1000.)
MD25=(data["boxLength"]==2500.)
MD40=(data["boxLength"]==4000.)

log_vmax_04 = (n.log10(data["log_"+qty+"_min"][ok & MD04])+n.log10(data["log_"+qty+"_max"][ok & MD04]))/2.
error_04 = data['std90_pc_cen'][ok & MD04]

log_vmax_10 = (n.log10(data["log_"+qty+"_min"][ok & MD10])+n.log10(data["log_"+qty+"_max"][ok & MD10]))/2.
error_10 = data['std90_pc_cen'][ok & MD10]

log_vmax_25 = (n.log10(data["log_"+qty+"_min"][ok & MD25]) +n.log10(data["log_"+qty+"_max"][ok & MD25]))/2.
error_25 = data['std90_pc_cen'][ok & MD25]

log_vmax_40 = (n.log10(data["log_"+qty+"_min"][ok & MD40])+n.log10(data["log_"+qty+"_max"][ok & MD40]))/2.
error_40 = data['std90_pc_cen'][ok & MD40]

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
p.plot(data["std90_pc_"+cos], data["dN_counts_"+cos]**(-0.5), 'ko', label='all', alpha=0.01)
p.plot(data["std90_pc_"+cos][MD04], data["dN_counts_"+cos][MD04]**(-0.5),marker='x',label="MD04",ls='')
p.plot(data["std90_pc_"+cos][MD10], data["dN_counts_"+cos][MD10]**(-0.5),marker='+',label="MD10",ls='')
p.plot(data["std90_pc_"+cos][MD25], data["dN_counts_"+cos][MD25]**(-0.5),marker='^',label="MD25",ls='')
p.plot(data["std90_pc_"+cos][MD40], data["dN_counts_"+cos][MD40]**(-0.5),marker='v',label="MD40",ls='')
xx = n.logspace(-4,0,20)
p.plot(xx, xx*3., ls='--', label='y=3x')
#p.axhline(Npmin**-0.5, c='r', ls='--', label='min counts cut')#r'$1/\sqrt{10^3}$')
p.axhline((10**6.87)**-0.5, c='k', ls='--', label='min vmax cut')#r'$1/\sqrt{10^{4.87}}$')
#p.xlim((2e-4,4e-1))
#p.ylim((2e-4,4e-1))
p.ylabel(r'$1/\sqrt{count} \; [\%]$')
p.xlabel(r'Jackknife  Resampling Error [%]')
p.yscale('log')
p.xscale('log')
gl = p.legend(loc=0,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig(join(dir,qty,"vmax-"+cos+"-jackknife-countsSqrt.png"))
p.clf()
from os.path import join
import numpy as n
import astropy.io.fits as fits
import os
import lib_functions_1pt as lib
import cPickle 

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
#Quantity studied
#=================
qty = "vmax"
# working directory
#=================
dir = join(os.environ['MULTIDARK_LIGHTCONE_DIR'], qty)

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
zmax = 2.3

cos = "cen"

tolerance = 0.06

#=================
# DATA
#=================
# loads summary file
data = fits.open( join(dir, "MD_"+qty+"_summary.fits"))[1].data
# redshift selection
zSel = lib.zSelection( data, zmin, zmax )
# mass selection
mSel = lib.mSelection(data, qty, limits_04, limits_10, limits_25,limits_40) 
# minimum number counts selection
nSel = lib.nSelection(data, NminCount, cos)
# altogether
ok = (zSel) & (mSel) & (nSel)


# x coordinates definition
#=================
log_vmax = (n.log10(data["log_"+qty+"_min"])+n.log10(data["log_"+qty+"_max"]))/2.
vmax = 10**log_vmax
x_data = log_vmax[ok]

# y coordinates
#=================
norm = (100)**3. /(cosmo.H(data["redshift"]).value)**6.
log_VF = n.log10( norm * vmax**3. * data["dNdVdlnM_"+cos])
log_VF_c = n.log10( norm * vmax**3. * data["dNdVdlnM_"+cos+"_c"])
y_data = log_VF[ok]

# z coordinate 
#=================
z_data = data["redshift"][ok]

# error on y position
#=================
error = data["dN_counts_"+cos]**(-0.5)
y_err = error[ok]
#error_04 = (data['std90_pc_cen'][MD04]*errorFactor+systError)/ n.log(10.) 


# FITTING THE TREND with REDSHIFT
#=======================
param_z0_file=open(join(dir, "vmax-"+cos+"-diff-function-z0-params.pkl"), 'r')
outCF = cPickle.load(param_z0_file)
param_z0_file.close()
A0, v0, a0, b0 = outCF[0]
print "----------------------------------------------------------"
print A0, v0, a0, b0
print "----------------------------------------------------------"
Az = lambda z, A1 : A0 + z*A1# + z**2 *A2
#vz = lambda z, v1, v2: v0 + z*v1 + v2*z**2.
vz = lambda z, v3 : v0 +v3*z**3.
az = lambda z, a3 : a0 + a3*z**1.
#bz = lambda z, b1, b2 : b0 +z*b1 + z**2.*b2
bz = lambda z, b1, b3 : b0 +z*b1 +b3*z**3.

#4.97 2.99 0.35 -0.17
#5.11 2.8 0.268 -0.15

ps = [ 0., 0., 0., 0., 0.]
"""
vf = lambda v, z, A1, v1, v2, a1, b1, b2 : n.log10( 10**Az(z, A1) * (10**v/10**vz(z, v1,v2))**(-bz(z, b1, b2)) * n.e**(- (10**v/10**vz(z, v1,v2))**(10**az(z, a1)) ) )
logFun = lambda v, z, ps : vf(v, z, ps[0], ps[1], ps[2], ps[3], ps[4], ps[5]) 
"""

vf = lambda v, z, A1,  v3, a3, b1, b3: n.log10( 10**Az(z, A1) * (10**v/10**vz(z,  v3))**(-bz(z, b1,b3)) * n.e**(- (10**v/10**vz(z,  v3))**(az(z,a3) ) ) )
logFun = lambda v, z, ps : vf(v, z, ps[0], ps[1], ps[2], ps[3], ps[4])#, ps[5], ps[6])

# chi2fun = lambda ps : n.sum( (logFun(x_data, ps) - y_data)**2. / (y_err)**2. )/(len(y_data) - len(ps))
chi2fun = lambda ps : n.sum( abs(logFun(x_data, z_data, ps) - y_data) / (y_err) )/(len(y_data) - len(ps))



res = minimize(chi2fun, ps, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 50000000000000})
pOpt = res.x
cov = res.direc
#chi2perpoint = lambda ps : (funG(lg_vmax, lg_1pz, ps) - lg_MF_c)**2. / (errorLog)**2. 
#chi2pp = chi2perpoint(pOpt)
#|print pOpt, cov
print "----------------------------------------------------------"
print n.round(pOpt,3)
print n.round(abs(cov.diagonal())**0.5,3)
print "----------------------------------------------------------"

outfile=open(join(dir, "vmax-"+cos+"-diff-function-z1-params.pkl"), 'w')
cPickle.dump(res, outfile)
outfile.close()

X = n.arange(n.min(x_data),n.max(x_data), 0.01)
Z = n.arange(zmin, zmax, 0.01)
x_model, z_model = n.meshgrid(X, Z)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])

#sc1=p.scatter(-n.hstack((x_model)), logFun(n.hstack((x_model)), n.hstack((z_model)), pOpt) , c=n.hstack((z_model)), s=2, marker='o',label="model", rasterized=True, vmin=zmin, vmax=zmax)

sc1=p.scatter(x_data, logFun(x_data, z_data, pOpt) , c=z_data, s=5, marker='o',label="model", rasterized=True, vmin=zmin, vmax = zmax)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'log$_{10}[V_{max}/(km \; s^{-1})]$')
p.ylabel(r'log$_{10} V^4 dn(V)/dlnV$') # log$_{10}[ n(>M)]')
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)
#p.yscale('log')
p.grid()
p.savefig(join(dir, "vmax-"+cos+"-differential-function-redshift-trend-model.png"))
p.clf()

"""

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(-x_data, y_data , c=z_data, s=5, marker='o',label="data", rasterized=True, vmin=zmin, vmax=zmax)
#p.errorbar(-x_data, y_data , yerr=y_err,fmt='none', rasterized=True)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'log$_{10}[V_{max}/(km \; s^{-1})]$')
p.ylabel(r'log$_{10} V^4 n(>V)$') # log$_{10}[ n(>M)]')
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)
p.xlim((-0.6, 0.4))
p.ylim((-3, 0.))
#p.yscale('log')
p.grid()
p.savefig(join(dir,qty,"vmax-"+cos+"-differential-function-redshift-trend-data.png"))
p.clf()

"""

MD_sel_fun=lambda name : (data["boxName"]==name)
MDnames= n.array(['MD_0.4Gpc', 'MD_1Gpc_new_rockS', 'MD_2.5Gpc','MD_4Gpc','MD_2.5GpcNW','MD_4GpcNW'])
MDsels=n.array([MD_sel_fun(name)[ok] for name in MDnames])

f_diff_fun = lambda MDs:  y_data[MDs] - logFun(x_data[MDs], z_data[MDs], pOpt)
f_diffs = n.array([f_diff_fun(MD) for MD in MDsels])

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
for index, fd in enumerate(f_diffs):
	inTol = (abs(10**fd-1)<tolerance)
	print index
	if len(fd)>0:
		p.errorbar(x_data[MDsels[index]], 10**fd, yerr = y_err[MDsels[index]] , rasterized=True, fmt='none', label=MDnames[index])
		print len(inTol.nonzero()[0]), len(fd), 100.*len(inTol.nonzero()[0])/ len(fd)

p.axhline(1.01,c='k',ls='--',label=r'syst $\pm1\%$')
p.axhline(0.99,c='k',ls='--')
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)

p.xlabel(r'log$_{10}[V_{max}/(km \; s^{-1})]$')
p.ylabel(r'data / model') 
#p.xlim((-0.6, 0.4))
#p.ylim((-0.9, 1.1))
#p.yscale('log')
p.grid()
p.savefig(join(dir,"vmax-"+cos+"-differential-function-redshift-trend-residual.png"))
p.clf()

"""
p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(x_data, f_diff , c=z_data, s=5, marker='o',label="model reiduals", rasterized=True, vmin=zmin, vmax = zmax)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)

p.xlabel(r'log$_{10}[V_{max}/(km \; s^{-1})]$')
p.ylabel(r'data - model') 
#p.xlim((-0.6, 0.4))
#p.ylim((-0.9, 1.1))
#p.yscale('log')
p.grid()
p.savefig(join(dir,qty,"vmax-"+cos+"-differential-function-redshift-trend-residual-z.png"))
p.clf()


tolerance = 0.06
in04 = (abs(10**f_diff_04-1)<tolerance)
print len(in04.nonzero()[0]), len(f_diff_04), 100.*len(in04.nonzero()[0])/ len(f_diff_04)
in10 = (abs(10**f_diff_10-1)<tolerance)
print len(in10.nonzero()[0]), len(f_diff_10), 100.*len(in10.nonzero()[0])/ len(f_diff_10)
in25 = (abs(10**f_diff_25-1)<tolerance)
print len(in25.nonzero()[0]), len(f_diff_25), 100.*len(in25.nonzero()[0])/ len(f_diff_25)
in40 = (abs(10**f_diff_40-1)<tolerance)
print len(in40.nonzero()[0]), len(f_diff_40), 100.*len(in40.nonzero()[0])/ len(f_diff_40)

tolerance = 0.08
inall = (abs(10**f_diff-1)<tolerance)
print len(inall.nonzero()[0]), len(f_diff), 100.*len(inall.nonzero()[0])/ len(f_diff)
"""
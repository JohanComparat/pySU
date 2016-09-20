from os.path import join
import numpy as n
import astropy.io.fits as fits
import os

import lib_functions_1pt as lib

import astropy.cosmology as co
cosmo = co.Planck13
import astropy.units as uu

#Quantity studied
qty = "mvir"
cos = 'cen'
# working directory
dir = join(os.environ['MULTIDARK_LIGHTCONE_DIR'], qty)
# loads summary file
data = fits.open( join(dir, "MD_"+qty+"_summary.fits"))[1].data

NminCount = 1000
Npmin = 1# 300
nolim = [0,1e17]
limits_04 = nolim # [Npmin*9.63 * 10**7, 5e12]
limits_10 = nolim # [Npmin*1.51 * 10**9., 5e13]
limits_25 = nolim # [Npmin*2.359 * 10**10., 5e14]
limits_40 = nolim # [Npmin* 9.6 * 10**10. , 5e15]
MPART = n.array([9.63 * 10**7, 1.51 * 10**9, 2.359 * 10**10, 9.6 * 10**10])

zmin = -0.01
zmax = 2.3

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
logsigM1 = n.log10(1./data['sigmaM'])#
log_mvir = (data["log_"+qty+"_min"]+data["log_"+qty+"_max"])/2.
mvir = 10**log_mvir

# y coordinates
rhom = cosmo.critical_density(data["redshift"]).to(uu.solMass/(uu.Mpc)**3.)/(cosmo.H(data["redshift"])/(100*uu.km/(uu.Mpc*uu.s)))**1.
log_MF = n.log10( mvir * data["dNdVdlnM_"+cos]/ rhom.value )
log_MF_c = n.log10(  data["dNdVdlnM_"+cos+"_c"])
log_f =  n.log10(mvir * data["dNdVdlnM_"+cos]/ rhom.value  / abs(data["dlnsigmaM1_o_dlnM"]))
log_f_c =  n.log10(mvir * data["dNdVdlnM_"+cos+"_c"]/ rhom.value  / abs(data["dlnsigmaM1_o_dlnM"]))

# NOW PLOTTING ALL THE DATA
lib.plot_mvir_function_data(log_mvir[ok], logsigM1[ok], log_MF[ok], log_MF_c[ok], data['redshift'][ok], zmin, zmax, cos = cos)

lib.plot_mvir_function_data_perBox(log_mvir, log_MF, MD04, MD10, MD25, MD25NW, MD40, MD40NW, cos=cos)

# ERROR PLOT: JK vs. POISSON
x = data["std90_pc_"+cos] 
y = data["dN_counts_"+cos]**(-0.5)
lib.plot_vmax_function_jackknife_poisson_error(x, y, MD04, MD10, MD25, MD25NW, MD40, MD40NW, cos = cos)


cos = 'sat'


nSelSat = lib.nSelection(data, NminCount, cos )
# altogether
ok = (zSel) & (mSel) & (nSelSat)
# selection per box :
MD04=(data["boxName"]=='MD_0.4Gpc')
MD10=(data["boxName"]=='MD_1Gpc_new_rockS')
MD25=(data["boxName"]=='MD_2.5Gpc')
MD40=(data["boxName"]=='MD_4Gpc')
MD25NW=(data["boxName"]=='MD_2.5GpcNW')
MD40NW=(data["boxName"]=='MD_4GpcNW')

# x coordinates definition
logsigM1 = n.log10(1./data['sigmaM'])#
log_mvir = (data["log_"+qty+"_min"]+data["log_"+qty+"_max"])/2.
mvir = 10**log_mvir

# y coordinates
rhom = cosmo.critical_density(data["redshift"]).to(uu.solMass/(uu.Mpc)**3.)/(cosmo.H(data["redshift"])/(100*uu.km/(uu.Mpc*uu.s)))**1.
log_MF = n.log10( mvir * data["dNdVdlnM_"+cos]/ rhom.value )
log_MF_c = n.log10(  data["dNdVdlnM_"+cos+"_c"])
log_f =  n.log10(mvir * data["dNdVdlnM_"+cos]/ rhom.value  / abs(data["dlnsigmaM1_o_dlnM"]))
log_f_c =  n.log10(mvir * data["dNdVdlnM_"+cos+"_c"]/ rhom.value  / abs(data["dlnsigmaM1_o_dlnM"]))

# NOW PLOTTING ALL THE DATA
lib.plot_mvir_function_data(log_mvir[ok], logsigM1[ok], log_MF[ok], log_MF_c[ok], data['redshift'][ok], zmin, zmax, cos = cos)

lib.plot_mvir_function_data_perBox(log_mvir, log_MF, MD04, MD10, MD25, MD25NW, MD40, MD40NW, cos=cos)

# ERROR PLOT: JK vs. POISSON
x = data["std90_pc_"+cos] 
y = data["dN_counts_"+cos]**(-0.5)
lib.plot_vmax_function_jackknife_poisson_error(x, y, MD04, MD10, MD25, MD25NW, MD40, MD40NW, cos = cos)

sys.exit()
"""
# PLOTTING THE ERROR PER BOX
lib.plot_vmax_function_data_error(log_vmax[ok & MD04], data['std90_pc_'+cos][ok & MD04], data["redshift"][ok & MD04], label='MD04', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data04-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD10], data['std90_pc_'+cos][ok & MD10], data["redshift"][ok & MD10], label='MD10', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data10-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD25], data['std90_pc_'+cos][ok & MD25], data["redshift"][ok & MD25], label='MD25', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data25-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD25NW], data['std90_pc_'+cos][ok & MD25NW], data["redshift"][ok & MD25NW], label='MD25NW', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data25NW-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD40], data['std90_pc_'+cos][ok & MD40], data["redshift"][ok & MD40], label='MD40', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data40-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD40NW], data['std90_pc_'+cos][ok & MD40NW], data["redshift"][ok & MD40NW], label='MD40NW', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data40NW-uncertainty.png")
"""
cos = 'sat'

# minimum number counts selection
nSelSat = lib.nSelection(data, NminCount, cos )
# altogether
ok = (zSel) & (mSel) & (nSelSat)

# y coordinates
log_VF = n.log10( norm * vmax**3. * data["dNdVdlnM_"+cos])
log_VF_c = n.log10( norm * vmax**3. * data["dNdVdlnM_"+cos+"_c"])
#print n.min(log_VF), n.max(log_VF)

# NOW PLOTTING ALL THE DATA
lib.plot_vmax_function_data(log_vmax, log_VF, log_VF_c, data["redshift"], zmin = -0.01, zmax = 2.3, cos=cos)

lib.plot_vmax_function_data_perBox(log_vmax, log_VF, log_VF_c, MD04, MD10, MD25, MD25NW, MD40, MD40NW, cos=cos)

"""
# PLOTTING THE ERROR PER BOX
lib.plot_vmax_function_data_error(log_vmax[ok & MD04], data['std90_pc_'+cos][ok & MD04], data["redshift"][ok & MD04], label='MD04', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data04-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD10], data['std90_pc_'+cos][ok & MD10], data["redshift"][ok & MD10], label='MD10', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data10-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD25], data['std90_pc_'+cos][ok & MD25], data["redshift"][ok & MD25], label='MD25', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data25-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD25NW], data['std90_pc_'+cos][ok & MD25NW], data["redshift"][ok & MD25NW], label='MD25NW', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data25NW-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD40], data['std90_pc_'+cos][ok & MD40], data["redshift"][ok & MD40], label='MD40', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data40-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD40NW], data['std90_pc_'+cos][ok & MD40NW], data["redshift"][ok & MD40NW], label='MD40NW', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data40NW-uncertainty.png")
"""
# ERROR PLOT: JK vs. POISSON
x = data["std90_pc_"+cos] 
y = data["dN_counts_"+cos]**(-0.5)
lib.plot_vmax_function_jackknife_poisson_error(x, y, MD04, MD10, MD25, MD25NW, MD40, MD40NW, cos = "cen")

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

from scipy.interpolate import interp1d
from scipy.misc import derivative

sigma = n.arange(0.05,10,0.05)
f = lambda sigma, A, a, b, c : A*( (sigma/b)**(-a) + 1 )*n.e**(-c/sigma**2.)
ft08Y = f(sigma, 0.224, 1.67, 1.80, 1.48) 
ft08X = n.log10(1./sigma)


dir='..'
dir_04 = join(dir,"MD_0.4Gpc")
dir_10 = join(dir,"MD_1Gpc")
dir_25 = join(dir,"MD_2.5Gpc")
dir_40 = join(dir,"MD_4Gpc")
dir_25N = join(dir,"MD_2.5GpcNW")
dir_40N = join(dir,"MD_4GpcNW")

data = fits.open( join("..", "mvir", "MD_mvir_summary.fits") )[1].data

boxes = set(data['boxName'])
mk = {"MD_0.4Gpc": '1', "MD_1Gpc": 'x' ,"MD_2.5Gpc": '2',"MD_4Gpc": '3' ,"MD_2.5GpcNW": '4', "MD_4GpcNW": '+'}

NminCount = 1 #10
Npmin = 1# 300
nolim = [0,1e17]
limits_04 = nolim # [Npmin*9.63 * 10**7, 5e12]
limits_10 = nolim # [Npmin*1.51 * 10**9., 5e13]
limits_25 = nolim # [Npmin*2.359 * 10**10., 5e14]
limits_40 = nolim # [Npmin* 9.6 * 10**10. , 5e15]
MPART = n.array([9.63 * 10**7, 1.51 * 10**9, 2.359 * 10**10, 9.6 * 10**10])
names = n.array(["SMD", "MDPL", "BigMD", "HMD", "BigMDNW", "HMDNW"])

zmin = -0.01
zmax = 0.001

def plotData(qty , cos , zmin = -0.01, zmax = 2.3):
	"""
	Plots the data to be used in the fits later in the analysis.
	"""
	# redshift selection
	zSel = (data["redshift"]>zmin)&(data["redshift"]<zmax)
	# mass selection
	mSel = ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>n.log10(limits_04[0])) &(data["log_"+qty+"_max"]<n.log10(limits_04[1]))) | ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>n.log10(limits_10[0])) &(data["log_"+qty+"_max"]<n.log10(limits_10[1]))) |  ((data["boxLength"]==2500.)&(data["log_"+qty+"_min"]>n.log10(limits_25[0])) &(data["log_"+qty+"_max"]<n.log10(limits_25[1]))) |  ((data["boxLength"]==4000.)&(data["log_"+qty+"_min"]>n.log10(limits_40[0]))&(data["log_"+qty+"_max"]<n.log10(limits_40[1]))) 
	#if  cos == "sat":
	#mSel =  ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>n.log10(limits_04[0])) &(data["log_"+qty+"_max"]<n.log10(limits_04[1]))) |  ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>n.log10(limits_10[0])) &(data["log_"+qty+"_max"]<n.log10(limits_10[1])))
	# minimum number counts selection
	nSel = (data['dN_counts_'+cos]>NminCount)
	# altogether
	ok = (zSel) & (mSel) & (nSel)
	# x coordinates
	logsigM1 = n.log10(1./data['sigmaM'])#
	print n.min(logsigM1), n.max(logsigM1)
	log_mvir = (data["log_"+qty+"_min"]+data["log_"+qty+"_max"])/2.
	mvir = 10**log_mvir
	# mean density array normalization
	rhom = cosmo.critical_density(data["redshift"]).to(uu.solMass/(uu.Mpc)**3.)/(cosmo.H(data["redshift"])/(100*uu.km/(uu.Mpc*uu.s)))**1.
	# y coordinates
	log_MF = n.log10( mvir * data["dNdVdlnM_"+cos]/ rhom.value )
	log_MF_c = n.log10(  data["dNdVdlnM_"+cos+"_c"])
	log_f =  n.log10(mvir * data["dNdVdlnM_"+cos]/ rhom.value  / abs(data["dlnsigmaM1_o_dlnM"]))
	log_f_c =  n.log10(mvir * data["dNdVdlnM_"+cos+"_c"]/ rhom.value  / abs(data["dlnsigmaM1_o_dlnM"]))
	
	# now the plots
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	for box in boxes:
		bb = (data["boxName"]==box)
		#print box, mk[box]
		sc1=p.scatter(logsigM1[bb], log_MF[bb], c=data["redshift"][bb], s=20, marker=mk[box],label=box, rasterized=True, vmin=zmin, vmax = zmax, edgecolors='face')
	
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'$ln(\sigma^{-1})$')
	p.ylabel(r'log$_{10} (M^2/\rho_m) dn(M)/dM$') 
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	p.xlim((-0.7,0.6))
	p.ylim((-4.5,-2))
	#p.yscale('log')
	p.grid()
	p.savefig(join(dir,qty,qty+"-"+cos+"-differential-function-data-xSigma.png"))
	p.clf()

	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	for box in boxes:
		bb = (data["boxName"]==box)
		#print box, mk[box]
		sc1=p.scatter(log_mvir[bb], log_MF[bb], c=data["redshift"][bb], s=20, marker=mk[box],label=box, rasterized=True, vmin=zmin, vmax = zmax, edgecolors='face')
	
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
	p.ylabel(r'log$_{10} [(M^2/\rho_m) dn(M)/dM]$') # log$_{10}[ n(>M)]')
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	#p.ylim((-8,1))
	p.xlim((9.5,16))
	p.ylim((-4.5,-2))
	p.grid()
	p.savefig(join(dir,qty,qty+"-"+cos+"-differential-function-data.png"))
	p.clf()
	"""
	MD04=(data["boxLength"]==400.)
	MD10=(data["boxLength"]==1000.)
	MD25=(data["boxLength"]==2500.)
	MD40=(data["boxLength"]==4000.)
	log_mvir_04 = (data["log_"+qty+"_min"][ok & MD04]+data["log_"+qty+"_max"][ok & MD04])/2.
	error_04 = data['std90_pc_cen'][ok & MD04]
	
	log_mvir_10 = (data["log_"+qty+"_min"][ok & MD10]+data["log_"+qty+"_max"][ok & MD10])/2.
	error_10 = data['std90_pc_cen'][ok & MD10]
	
	log_mvir_25 = (data["log_"+qty+"_min"][ok & MD25]+data["log_"+qty+"_max"][ok & MD25])/2.
	error_25 = data['std90_pc_cen'][ok & MD25]
	
	log_mvir_40 = (data["log_"+qty+"_min"][ok & MD40]+data["log_"+qty+"_max"][ok & MD40])/2.
	error_40 = data['std90_pc_cen'][ok & MD40]
	
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	sc1=p.scatter(log_mvir_04, 100*error_04, c=data["redshift"][ok & MD04], s=5, marker='o',label="MD04", rasterized=True, vmin=zmin, vmax = zmax)
	sc1.set_edgecolor('face')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
	p.ylabel(r'Jackknife resampling relative error [%]') # log$_{10}[ n(>M)]')
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	p.ylim((2e-2,30))
	p.xlim((9.5,16))
	p.yscale('log')
	p.grid()
	p.savefig(join(dir,qty,qty+"-"+cos+"-cumulative-function-data04-uncertainty.png"))
	p.clf()
	
	
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	sc1=p.scatter(log_mvir_10, 100*error_10, c=data["redshift"][ok & MD10], s=5, marker='o',label="MD10", rasterized=True, vmin=zmin, vmax = zmax)
	sc1.set_edgecolor('face')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
	p.ylabel(r'Jackknife resampling relative error [%]') # log$_{10}[ n(>M)]')
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	p.ylim((2e-2,30))
	p.xlim((9.5,16))
	p.yscale('log')
	p.grid()
	p.savefig(join(dir, qty, qty+"-"+cos+"-cumulative-function-data10-uncertainty.png"))
	p.clf()
	
	
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	sc1=p.scatter(log_mvir_25, 100*error_25, c=data["redshift"][ok & MD25], s=5, marker='o',label="MD25", rasterized=True, vmin=zmin, vmax = zmax)
	sc1.set_edgecolor('face')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
	p.ylabel(r'Jackknife resampling relative error [%]') # log$_{10}[ n(>M)]')
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	p.ylim((2e-2,30))
	p.xlim((9.5,16))
	p.yscale('log')
	p.grid()
	p.savefig(join(dir,qty,qty+"-"+cos+"-cumulative-function-data25-uncertainty.png"))
	p.clf()
	
	
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	sc1=p.scatter(log_mvir_40, 100*error_40, c=data["redshift"][ok & MD40], s=5, marker='o',label="MD40", rasterized=True, vmin=zmin, vmax = zmax)
	sc1.set_edgecolor('face')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
	p.ylabel(r'Jackknife resampling relative error [%]') # log$_{10}[ n(>M)]')
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	p.ylim((2e-2,30))
	p.xlim((9.5,16))
	p.yscale('log')
	p.grid()
	p.savefig(join(dir,qty,qty+"-"+cos+"-cumulative-function-data40-uncertainty.png"))
	p.clf()
	
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	sc1=p.scatter(log_mvir, log_MF_c, c=data["redshift"], s=5, marker='o',label="data", rasterized=True, vmin=zmin, vmax = zmax)
	sc1.set_edgecolor('face')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
	p.ylabel(r'log$_{10} n(>M)$') # log$_{10}[ n(>M)]')
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	#p.ylim((-8,1))
	p.xlim((9.5,16))
	#p.yscale('log')
	p.grid()
	p.savefig(join(dir,qty,qty+"-"+cos+"-cumulative-function-data.png"))
	p.clf()
	
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	sc1=p.scatter(logsigM1, log_f_c, c=data["redshift"], s=5, marker='o',label="data", rasterized=True, vmin=zmin, vmax = zmax)
	sc1.set_edgecolor('face')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'$log(\sigma^{-1})$')
	p.ylabel(r'log$_{10} n(>M)$') # log$_{10}[ n(>M)]')
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	#p.ylim((-8,1))
	#p.xlim((9.5,16))
	#p.xlim((-0.7,0.6))
	#p.ylim((-5.5,2))
	#p.yscale('log')
	p.grid()
	p.savefig(join(dir,qty,qty+"-"+cos+"-cumulative-function-data-xSigma.png"))
	p.clf()
	"""
plotData(qty = 'mvir', cos = "sat")

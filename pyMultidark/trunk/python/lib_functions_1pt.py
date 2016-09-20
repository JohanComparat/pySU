import glob
import sys
import astropy.io.fits as fits
import os
from os.path import join
import cPickle

import numpy as n
from scipy.interpolate import interp1d
from scipy.misc import derivative

import matplotlib
matplotlib.use('pdf')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p

# MULTIDARK TABLE GENERIC FUNCTIONS
vf = lambda v, A, v0, alpha, beta : n.log10( 10**A * (10**v/10**v0)**(-beta) * n.e**(- (10**v/10**v0)**(alpha) ) )

mSelection = lambda data, qty, limits_04, limits_10, limits_25, limits_40 : ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>limits_04[0]) &(data["log_"+qty+"_max"]<limits_04[1])) | ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>limits_10[0]) &(data["log_"+qty+"_max"]<limits_10[1])) |  ((data["boxLength"]==2500.)&(data["log_"+qty+"_min"]>limits_25[0]) &(data["log_"+qty+"_max"]<limits_25[1])) |  ((data["boxLength"]==4000.)&(data["log_"+qty+"_min"]>limits_40[0])&(data["log_"+qty+"_max"]<limits_40[1])) 

zSelection = lambda data, zmin, zmax : (data["redshift"]>zmin)&(data["redshift"]<zmax)

nSelection = lambda data, NminCount, cos : (data['dN_counts_'+cos]>NminCount)

# VMAX 1point FUNCTION 
def plot_vmax_function_jackknife_poisson_error(x, y, MD04, MD10, MD25, MD25NW, MD40, MD40NW, cos = "cen", dir=join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'vmax')):
	"""
	:param x: x coordinates
	:param y: y coordinates
	:param cos: centra or satelitte. Default: "cen"
	:param dir: working directory. Default: join( os.environ['MULTIDARK_LIGHTCONE_DIR'], qty), :param qty: quantity studied. Default: 'vmax'
	"""
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	p.plot(x, y, 'ko', label='all', alpha=0.01)
	p.plot(x[MD04], y[MD04]**(-0.5),marker='x',label="MD04",ls='')
	p.plot(x[MD10], y[MD10]**(-0.5),marker='+',label="MD10",ls='')
	p.plot(x[MD25], y[MD25]**(-0.5),marker='^',label="MD25",ls='')
	p.plot(x[MD40], y[MD40]**(-0.5),marker='v',label="MD40",ls='')
	p.plot(x[MD25NW], y[MD25NW]**(-0.5),marker='^',label="MD25NW",ls='')
	p.plot(x[MD40NW], y[MD40NW]**(-0.5),marker='v',label="MD40NW",ls='')
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
	p.savefig(join(dir,"vmax-"+cos+"-jackknife-countsSqrt.png"))
	p.clf()

def plot_vmax_function_data_error(log_vmax, error, redshift, label, zmin, zmax, cos = "cen", figName="vmax-cen-data04-uncertainty.png", dir=join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'vmax')):
	"""
	:param log_vmax: x coordinates
	:param error: y coordinates
	:param redshift: color coordinate
	:param label: label in the caption
	:param zmin: minimum redshift
	:param zmax: maximum redshift
	:param cos: centra or satelitte. Default: "cen"
	:param figName: string to be added to the figure name. Default:=""
	:param dir: working directory. Default: join( os.environ['MULTIDARK_LIGHTCONE_DIR'], qty), :param qty: quantity studied. Default: 'vmax'
	"""
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	sc1=p.scatter(log_vmax, 100*error, c=redshift, s=5, marker='o',label=label, rasterized=True, vmin=zmin, vmax = zmax)
	sc1.set_edgecolor('face')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[V_{max}/(km \; s^{-1})]$')
	p.ylabel(r'JK relative error [%]') 
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	p.ylim((2e-2,30))
	#p.xlim((1.5, 3.5))
	p.yscale('log')
	p.grid()
	p.savefig(join(dir,figName))
	p.clf()

def plot_vmax_function_data(log_vmax, log_VF, log_VF_c, redshift, zmin, zmax, cos = "cen", figName="", dir=join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'vmax')):
	"""
	:param log_vmax: x coordinates
	:param log_VF: y coordinates
	:param redshift: color coordinate
	:param zmin: minimum redshift
	:param zmax: maximum redshift
	:param cos: centra or satelitte. Default: "cen"
	:param figName: string to be added to the figure name. Default:=""
	:param dir: working directory. Default: join( os.environ['MULTIDARK_LIGHTCONE_DIR'], qty), :param qty: quantity studied. Default: 'vmax'
	"""
	# now the plots
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	sc1=p.scatter(log_vmax, log_VF, c=redshift, s=5, marker='o',label="MD "+cos+" data", rasterized=True, vmin=zmin, vmax = zmax)
	sc1.set_edgecolor('face')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[V_{max}/(km \; s^{-1})]$')
	p.ylabel(r'log$_{10} [(V^3/H^3(z)\; dn(V)/dlnV]$') # log$_{10}[ n(>M)]')
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	#p.ylim((-8,1))
	#p.xlim((1.5, 3.5))
	#p.ylim((-3.5,-1))
	p.grid()
	p.savefig(join(dir,"vmax-"+figName+cos+"-differential-function-data.png"))
	p.clf()
	
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	sc1=p.scatter(log_vmax, log_VF_c, c=data["redshift"][ok], s=5, marker='o',label="data", rasterized=True, vmin=zmin, vmax = zmax)
	sc1.set_edgecolor('face')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[V_{max}/(km \; s^{-1})]$')
	p.ylabel(r'log$_{10} [V^3/H^3(z)\; n(>V)]$') # log$_{10}[ n(>M)]')
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	#p.ylim((-8,1))
	#p.xlim((1.5, 3.5))
	#p.yscale('log')
	p.grid()
	p.savefig(join(dir,"vmax-"+figName+cos+"-cumulative-function-data.png"))
	p.clf()

def fit_vmax_function_z0(data, x_data, y_data , y_err, p0, 	tolerance = 0.03, cos = "cen", mode = "curve_fit", dir=join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'vmax')):
	"""
	Fits a function to the vmax data
	:param data: data table of the selected points for the fit
	:param x_data: x coordinate
	:param y_data: y coordinate
	:param y_err: error
	:param p0: first guess
	:param tolerance: percentage error tolerance to compute how many points are outside of the fit
	:param cos: central or satelitte
	:param mode: fitting mode, "curve_fit" or "minimize"
	:param dir: working dir
	:param qty: vmax here
	:return: result of the fit: best parameter array and covariance matrix
	produces a plot of the residuals
	"""
	if mode == "curve_fit":
		print "mode: curve_fit"
		pOpt, pCov=curve_fit(vf, x_data, y_data, p0, y_err)#, bounds=boundaries)
		print "best params=",pOpt[0], pOpt[1], pOpt[2], pOpt[3]
		print "err=",pCov[0][0]**0.5, pCov[1][1]**0.5, pCov[2][2]**0.5, pCov[3][3]**0.5
		
	if mode == "minimize":
		print "mode: minimize"
		chi2fun = lambda ps : n.sum( (vf(x_data, ps) - y_data)**2. / (y_err)**2. )/(len(y_data) - len(ps))
		res = minimize(chi2fun, p0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
		pOpt = res.x
		pCov = res.direc
		print "best params=",pOpt[0], pOpt[1], pOpt[2], pOpt[3]
		print "err=",pCov[0][0]**0.5, pCov[1][1]**0.5, pCov[2][2]**0.5, pCov[3][3]**0.5
		
	x_model = n.arange(n.min(x_data),n.max(x_data),0.005)
	y_model = vf(x_model, outCF[0][0], outCF[0][1], outCF[0][2], outCF[0][3])
	n.savetxt(join(dir,"vmax-"+cos+"-differential-function-z0-model-pts.txt"),n.transpose([x_model, y_model]) )
	outfile=open(join(dir,"vmax-"+cos+"-diff-function-z0-params.pkl"), 'w')
	cPickle.dump([pOpt, pCov], outfile)
	outfile.close()
			
	f_diff =  y_data - vf(x_data, pOpt[0], pOpt[1], pOpt[2], pOpt[3])
	
	MD04=(data["boxName"]=='MD_0.4Gpc')
	MD10=(data["boxName"]=='MD_1Gpc_new_rockS')
	MD25=(data["boxName"]=='MD_2.5Gpc')
	MD40=(data["boxName"]=='MD_4Gpc')
	MD25NW=(data["boxName"]=='MD_2.5GpcNW')
	MD40NW=(data["boxName"]=='MD_4GpcNW')

	f_diff_04 =  y_data[MD04] - vf(x_data[MD04], pOpt[0], pOpt[1], pOpt[2], pOpt[3])
	f_diff_10 =  y_data[MD10] - vf(x_data[MD10], pOpt[0], pOpt[1], pOpt[2], pOpt[3])
	f_diff_25 =  y_data[MD25] - vf(x_data[MD25], pOpt[0], pOpt[1], pOpt[2], pOpt[3])
	f_diff_25NW =  y_data[MD25NW] - vf(x_data[MD25NW], pOpt[0], pOpt[1], pOpt[2], pOpt[3])
	f_diff_40 =  y_data[MD40] - vf(x_data[MD40], pOpt[0], pOpt[1], pOpt[2], pOpt[3])
	f_diff_40NW =  y_data[MD40NW] - vf(x_data[MD40NW], pOpt[0], pOpt[1], pOpt[2], pOpt[3])
	f_diffs = [f_diff_04, f_diff_10,f_diff_25, f_diff_25NW, f_diff_40, f_diff_40NW]
	
	print "================================"
	for fd in f_diffs:
		in04 = (abs(10**fd-1)<tolerance)
		print len(in04.nonzero()[0]), len(fd), 100.*len(in04.nonzero()[0])/ len(fd)
	
	# now the plots
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	p.errorbar(x_data_04, 10**f_diff_04, yerr = error_04 , rasterized=True, fmt='none', label="MD04")
	p.errorbar(x_data_10, 10**f_diff_10, yerr = error_10 , rasterized=True, fmt='none', label="MD10")
	p.errorbar(x_data_25, 10**f_diff_25, yerr = error_25 , rasterized=True, fmt='none', label="MD25")
	p.errorbar(x_data_40, 10**f_diff_40, yerr = error_40 , rasterized=True, fmt='none', label="MD40")
	p.errorbar(x_data_25NW, 10**f_diff_25NW, yerr = error_25NW , rasterized=True, fmt='none', label="MD25")
	p.errorbar(x_data_40NW, 10**f_diff_40NW, yerr = error_40NW , rasterized=True, fmt='none', label="MD40")
	p.axhline(1.01,c='k',ls='--',label=r'syst $\pm1\%$')
	p.axhline(0.99,c='k',ls='--')
	p.xlabel(r'$log(V_{max})$')
	p.ylabel(r'data/model') 
	gl = p.legend(loc=0,fontsize=10)
	gl.set_frame_on(False)
	#p.xlim((-0.7,0.6))
	#p.ylim((-0.05,0.05))
	#p.yscale('log')
	p.grid()
	p.savefig(join(dir,"vmax-"+cos+"-differential-function-fit-residual-log.png"))
	p.clf()
	return pOpt, pCov

# MULTIDARK DATA OUTPUT HANDLING
def getStat(file,volume,unitVolume):
	"""
	From the pickle file output by the Multidark class, we output the number counts (differential and cumulative) per unit volume per mass bin.
	:param file: filename
	:param volume: total volume of the box
	:param unitVolume: sub volume used in the jackknife
	:return: Number counts, cumulative number counts, count density, cumulative count density, jackknife mean, jackknife std, cumulative jackknife mean, cumulative jackknife std
	"""
	# print file
	data=cPickle.load(open(file,'r'))
	data_c = n.array([n.array([ n.sum(el[ii:]) for ii in range(len(el)) ]) for el in data])
	Ncounts = data.sum(axis=0) 
	Ncounts_c = data_c.sum(axis=0) # n.array([ n.sum(Ncounts[ii:]) for ii in range(len(Ncounts)) ])
	Nall = Ncounts / volume
	Nall_c = Ncounts_c / volume
	index=n.arange(int(data.shape[0]))
	n.random.shuffle( index )
	Ntotal = int(data.shape[0])
	# discard 100
	def get_mean_std( pcDiscard = 0.1):
		"""
		retrieves the mean and std from the jackknife
		:param pcDiscard:percentage to discard for the jackknife
		:return: mean, std, cumulative mean90 and cumulative std
		"""
		Ndiscard = Ntotal * pcDiscard
		resamp = n.arange(0,Ntotal+1, Ndiscard)
		N90 = n.array([n.sum(data[n.delete(n.arange(Ntotal), index[resamp[i]:resamp[i+1]])], axis=0) for i in range(len(resamp)-1)]) / (unitVolume*(Ntotal - Ndiscard) )
		mean90 = n.mean(N90, axis=0)
		std90 = n.std(N90, axis=0) / mean90
		N90_c = n.array([n.sum(data_c[n.delete(n.arange(Ntotal), index[resamp[i]:resamp[i+1]])], axis=0) for i in range(len(resamp)-1)]) / (unitVolume*(Ntotal - Ndiscard) )
		mean90_c = n.mean(N90_c, axis=0)
		std90_c = n.std(N90_c, axis=0) / mean90_c
		return mean90, std90, mean90_c, std90_c

	mean90, std90, mean90_c, std90_c = get_mean_std(0.1)
	#mean99, std99, mean99_c, std99_c = getMS(0.01)
	sel = Nall>1/volume
	# print std99[sel]/std90[sel]
	# print mean90[sel]/mean99[sel]
	# print Nall[sel]/mean99[sel]
	# print Nall[sel]/mean90[sel]
	return Ncounts, Ncounts_c, Nall, Nall_c, mean90, std90, mean90_c, std90_c

def convert_pkl_mass(fileC, fileS, binFile, zList_files,z0, z0short, qty='mvir'):
	"""
	:param qty: one point function variable.Default: mvir.
	:param fileC: file with the central halo statistics
	:param fileS: file with the satelitte halo statistics
	:param binFile: file with the bins
	:param zList_files: list of file with linking snapshot number and redshift
	:param z0: redshift - number relation to link to the files containing the linear theory (sigma M relation, P(k) ...)
	:param z0short: same as z0 but shorter
	:return: a fits table containing the one point function histograms
	"""
	boxName = fileC.split('/')[6]
	boxZN = float(fileC.split('/')[-1].split('_')[1])
	bins = n.loadtxt(binFile)
	dX = ( 10**bins[1:]  - 10**bins[:-1] ) #* n.log(10)
	dlnbin = dX / (10**(( bins[1:]  + bins[:-1] )/2.))
	print boxName
	if boxName=='MD_0.4Gpc' :
		boxLength = 400.
		boxRedshift = 1./boxZN - 1.
		logmp = n.log10(9.63 * 10**7)
		
	if boxName=='MD_1Gpc_new_rockS' :
		boxLength = 1000.
		boxRedshift = 1./boxZN - 1.
		logmp = n.log10(1.51 * 10**9)

	if boxName=='MD_2.5Gpc' :
		boxLength = 2500.
		nSN, aSN = n.loadtxt(zList_files[2], unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
		conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))
		boxRedshift =  conversion[boxZN] 
		logmp = n.log10(2.359 * 10**10)

	if boxName=='MD_4Gpc' :
		boxLength = 4000.
		nSN, aSN = n.loadtxt(zList_files[3], unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
		conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))
		boxRedshift =  conversion[boxZN] 
		logmp = n.log10(9.6 * 10**10)

	if boxName=='MD_2.5GpcNW' :
		boxLength = 2500.
		nSN, aSN = n.loadtxt(zList_files[4], unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
		conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))
		boxRedshift =  conversion[boxZN] 
		logmp = n.log10(2.359 * 10**10)

	if boxName=='MD_4GpcNW' :
		boxLength = 4000.
		nSN, aSN = n.loadtxt(zList_files[5], unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
		conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))
		boxRedshift =  conversion[boxZN] 
		logmp = n.log10(9.6 * 10**10)

	index = int(n.argwhere( abs(z0-n.round(boxRedshift, 6))<0.00001)[0] )
	msigmaFile=join(os.environ['PYSU_MD_DIR'], "data", "PK_DM_CLASS", "hmf_highz_medz_lowz_planck", "mVector_z_"+str(z0short[index])+".txt")
	# print boxRedshift
	# print msigmaFile
	DATA = n.loadtxt(msigmaFile,unpack=True)
	# [1] m:            [M_sun/h] 
	# [2] sigma 
	# [3] ln(1/sigma) 
	# [4] n_eff 
	# [5] f(sigma) 
	# [6] dn/dm:        [h^4/(Mpc^3*M_sun)] 
	# [7] dn/dlnm:      [h^3/Mpc^3] 
	# [8] dn/dlog10m:   [h^3/Mpc^3] 
	# [9] n(>m):        [h^3/Mpc^3] 
	# [11] rho(>m):     [M_sun*h^2/Mpc^3] 
	# [11] rho(<m):     [M_sun*h^2/Mpc^3] 
	# [12] Lbox(N=1):   [Mpc/h]
	# print index, z0[ int(index) ], boxRedshift
	M=DATA[0]
	sigma = DATA[1]
	m2sigma = interp1d(M, sigma)
	toderive = interp1d(n.log(M), DATA[2])
	
	# print n.min(M),n.max(M), n.min(bins), n.max(bins)
	ok = (bins[1:] > logmp+1.0)&(bins[1:]<15.2)
	sig = m2sigma( 10**((bins[:-1][ok]+bins[1:][ok])/2.) )
	# print '================='
	# print n.min(n.log(10**((bins[:-1][ok]+bins[1:][ok])/2.))),n.max(n.log(10**((bins[:-1][ok]+bins[1:][ok])/2.))), n.min(n.log(M)), n.max(n.log(M))
	#selX = (M > 5.*n.min(M)) &(M<n.max(M)/5.)
	dlnsigdm = derivative(toderive, n.log(10**((bins[:-1][ok]+bins[1:][ok])/2.)))

	col000 = fits.Column( name="boxName",format="14A", array= n.array([boxName for i in range(len(bins[:-1]))]))
	col0 = fits.Column( name="sigmaM",format="D", array= sig )
	col00 = fits.Column( name="nuM",format="D", array= 1.686/sig )
	col01 = fits.Column( name="dlnsigmaM1_o_dlnM",format="D", array= dlnsigdm )
	col1 = fits.Column( name="log_"+qty+"_min",format="D", array= bins[:-1][ok] )
	col2 = fits.Column( name="log_"+qty+"_max",format="D", array= bins[1:][ok] )
	col3 = fits.Column( name="redshift",format="D", array= boxRedshift * n.ones_like(bins[:-1][ok]) )
	col4 = fits.Column( name="boxLength",format="D", array= boxLength * n.ones_like(bins[:-1][ok]) )
	col4_2 = fits.Column( name="Mpart",format="D", array=  logmp* n.ones_like(bins[:-1][ok]) )

	unitVolume =  (boxLength*0.10)**3.
	volume = (boxLength)**3.

	Ncounts, Ncounts_c, Nall, Nall_c, mean90, std90, mean90_c, std90_c = getStat(fileC,volume,unitVolume)

	col5 = fits.Column( name="dN_counts_cen",format="D", array= Ncounts[ok] )
	col6 = fits.Column( name="dN_counts_cen_c",format="D", array= Ncounts_c[ok])
	col7 = fits.Column( name="dNdV_cen",format="D", array= Nall[ok] )
	col8 = fits.Column( name="dNdV_cen_c",format="D", array= Nall_c[ok] )
	col9 = fits.Column( name="dNdVdlnM_cen",format="D", array= Nall[ok]/dlnbin[ok] )
	col10 = fits.Column( name="dNdVdlnM_cen_c",format="D", array= Nall_c[ok]/dlnbin[ok] )
	col11 = fits.Column( name="std90_pc_cen",format="D", array= std90[ok] )
	col12 = fits.Column( name="std90_pc_cen_c",format="D", array= std90_c[ok] )

	Ncounts, Ncounts_c, Nall, Nall_c, mean90, std90, mean90_c, std90_c = getStat(fileS,volume,unitVolume)

	col5_s = fits.Column( name="dN_counts_sat",format="D", array= Ncounts[ok] )
	col6_s = fits.Column( name="dN_counts_sat_c",format="D", array= Ncounts_c[ok])
	col7_s = fits.Column( name="dNdV_sat",format="D", array= Nall[ok] )
	col8_s = fits.Column( name="dNdV_sat_c",format="D", array= Nall_c[ok] )
	col9_s = fits.Column( name="dNdVdlnM_sat",format="D", array= Nall[ok]/dlnbin[ok] )
	col10_s = fits.Column( name="dNdVdlnM_sat_c",format="D", array= Nall_c[ok] / dlnbin[ok] )
	col11_s = fits.Column( name="std90_pc_sat",format="D", array= std90[ok] )
	col12_s = fits.Column( name="std90_pc_sat_c",format="D", array= std90_c[ok] )


	hdu2 = fits.BinTableHDU.from_columns([col000, col1, col2, col3, col4, col4_2, col5, col6, col7, col8, col9, col10, col11, col12, col5_s, col6_s, col7_s, col8_s, col9_s, col10_s, col11_s, col12_s, col0, col00, col01])
	
	writeName = join(os.environ['MULTIDARK_LIGHTCONE_DIR'], qty, "data", "MD_"+boxName+"_"+str(boxRedshift)+"_"+qty+".fits")
	os.system("rm -rf "+ writeName)
	hdu2.writeto(writeName)

def convert_pkl_velocity(fileC, fileS, binFile, zList_files, qty='vmax'):
	"""
	:param qty: one point function variable. Default vmax.
	:param fileC: file with the central halo statistics
	:param fileS: file with the satelitte halo statistics
	:param binFile: file with the bins
	:param zList_files: list of file with linking snapshot number and redshift
	:return: a fits table containing the one point function histograms
	"""
	print fileC.split('/')[6]
	boxName = fileC.split('/')[6]
	boxZN = float(fileC.split('/')[-1].split('_')[1])
	bins = n.loadtxt(binFile)
	#10**n.arange(0,3.5,0.01)
	dX = ( bins[1:]  - bins[:-1] ) #* n.log(10)
	dlnbin = dX / (( bins[1:]  + bins[:-1] )/2.)
	print boxName
	if boxName=='MD_0.4Gpc' :
		boxLength = 400.
		boxRedshift = 1./boxZN - 1.
		logmp = n.log10(9.63 * 10**7)
		
	if boxName=='MD_1Gpc_new_rockS' :
		boxLength = 1000.
		boxRedshift = 1./boxZN - 1.
		logmp = n.log10(1.51 * 10**9)

	if boxName=='MD_2.5Gpc' :
		boxLength = 2500.
		nSN, aSN = n.loadtxt(zList_files[2], unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
		conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))
		boxRedshift =  conversion[boxZN] 
		logmp = n.log10(2.359 * 10**10)

	if boxName=='MD_4Gpc' :
		boxLength = 4000.
		nSN, aSN = n.loadtxt(zList_files[3], unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
		conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))
		boxRedshift =  conversion[boxZN] 
		logmp = n.log10(9.6 * 10**10)

	if boxName=='MD_2.5GpcNW' :
		boxLength = 2500.
		nSN, aSN = n.loadtxt(zList_files[4], unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
		conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))
		boxRedshift =  conversion[boxZN] 
		logmp = n.log10(2.359 * 10**10)

	if boxName=='MD_4GpcNW' :
		boxLength = 4000.
		nSN, aSN = n.loadtxt(zList_files[5], unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
		conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))
		boxRedshift =  conversion[boxZN] 
		logmp = n.log10(9.6 * 10**10)

	col0 = fits.Column( name="boxName",format="14A", array= n.array([boxName for i in range(len(bins[:-1]))]))
	col1 = fits.Column( name="log_"+qty+"_min",format="D", array= bins[:-1] )
	col2 = fits.Column( name="log_"+qty+"_max",format="D", array= bins[1:] )
	col3 = fits.Column( name="redshift",format="D", array= boxRedshift * n.ones_like(bins[:-1]) )
	col4 = fits.Column( name="boxLength",format="D", array= boxLength * n.ones_like(bins[:-1]) )
	col4_2 = fits.Column( name="Mpart",format="D", array=  logmp* n.ones_like(bins[:-1]) )

	unitVolume =  (boxLength*0.10)**3.
	volume = (boxLength)**3.

	Ncounts, Ncounts_c, Nall, Nall_c, mean90, std90, mean90_c, std90_c = getStat(fileC,volume,unitVolume)

	col5 = fits.Column( name="dN_counts_cen",format="D", array= Ncounts )
	col6 = fits.Column( name="dN_counts_cen_c",format="D", array= Ncounts_c)
	col7 = fits.Column( name="dNdV_cen",format="D", array= Nall )
	col8 = fits.Column( name="dNdV_cen_c",format="D", array= Nall_c )
	col9 = fits.Column( name="dNdVdlnM_cen",format="D", array= Nall/dlnbin )
	col10 = fits.Column( name="dNdVdlnM_cen_c",format="D", array= Nall_c/dlnbin )
	col11 = fits.Column( name="std90_pc_cen",format="D", array= std90 )
	col12 = fits.Column( name="std90_pc_cen_c",format="D", array= std90_c )

	Ncounts, Ncounts_c, Nall, Nall_c, mean90, std90, mean90_c, std90_c = getStat(fileS,volume,unitVolume)

	col5_s = fits.Column( name="dN_counts_sat",format="D", array= Ncounts )
	col6_s = fits.Column( name="dN_counts_sat_c",format="D", array= Ncounts_c)
	col7_s = fits.Column( name="dNdV_sat",format="D", array= Nall )
	col8_s = fits.Column( name="dNdV_sat_c",format="D", array= Nall_c )
	col9_s = fits.Column( name="dNdVdlnM_sat",format="D", array= Nall/dlnbin )
	col10_s = fits.Column( name="dNdVdlnM_sat_c",format="D", array= Nall_c / dlnbin )
	col11_s = fits.Column( name="std90_pc_sat",format="D", array= std90 )
	col12_s = fits.Column( name="std90_pc_sat_c",format="D", array= std90_c )


	hdu2 = fits.BinTableHDU.from_columns([col0, col1, col2, col3, col4, col4_2, col5, col6, col7, col8, col9, col10, col11, col12, col5_s, col6_s, col7_s, col8_s, col9_s, col10_s, col11_s, col12_s])
	
	writeName = join(os.environ['MULTIDARK_LIGHTCONE_DIR'], qty, "data", "MD_"+boxName+"_"+str(boxRedshift)+"_"+qty+".fits")
	os.system("rm -rf "+ writeName)
	hdu2.writeto(writeName)


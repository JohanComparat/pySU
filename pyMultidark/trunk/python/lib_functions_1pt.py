import glob
import sys
from os.path import join
import numpy as n
import astropy.io.fits as fits
import os
from scipy.interpolate import interp1d
from scipy.misc import derivative
import cPickle

def getStat(file,volume,unitVolume):
	"""
	From the pickle file output by the Multidark class, we output the number counts (differential and cumulative) per unit volume per mass bin.
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
	def getMS( pcDiscard = 0.1):
		Ndiscard = Ntotal * pcDiscard
		resamp = n.arange(0,Ntotal+1, Ndiscard)
		N90 = n.array([n.sum(data[n.delete(n.arange(Ntotal), index[resamp[i]:resamp[i+1]])], axis=0) for i in range(len(resamp)-1)]) / (unitVolume*(Ntotal - Ndiscard) )
		mean90 = n.mean(N90, axis=0)
		std90 = n.std(N90, axis=0) / mean90
		N90_c = n.array([n.sum(data_c[n.delete(n.arange(Ntotal), index[resamp[i]:resamp[i+1]])], axis=0) for i in range(len(resamp)-1)]) / (unitVolume*(Ntotal - Ndiscard) )
		mean90_c = n.mean(N90_c, axis=0)
		std90_c = n.std(N90_c, axis=0) / mean90_c
		return mean90, std90, mean90_c, std90_c

	mean90, std90, mean90_c, std90_c = getMS(0.1)
	#mean99, std99, mean99_c, std99_c = getMS(0.01)
	sel = Nall>1/volume
	# print std99[sel]/std90[sel]
	# print mean90[sel]/mean99[sel]
	# print Nall[sel]/mean99[sel]
	# print Nall[sel]/mean90[sel]
	return Ncounts, Ncounts_c, Nall, Nall_c, mean90, std90, mean90_c, std90_c

def convert_pkl_mass(fileC, fileS, binFile, zList_files,z0, z0short, qty):
	"""returns a fits table containing the histograms
	"""
	boxName = fileC.split('/')[5]
	boxZN = float(fileC.split('/')[-1].split('_')[1])
	bins = n.loadtxt(binFile)
	dX = ( 10**bins[1:]  - 10**bins[:-1] ) #* n.log(10)
	dlnbin = dX / (10**(( bins[1:]  + bins[:-1] )/2.))

	if boxName=='MD_0.4Gpc' :
		boxLength = 400.
		boxRedshift = 1./boxZN - 1.
		logmp = n.log10(9.63 * 10**7)
		
	if boxName=='MD_1Gpc' :
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
	msigmaFile=join("..", "Pk_DM_CLASS", "hmf_highz_medz_lowz_planck", "mVector_z_"+str(z0short[index])+".txt")
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
	if os.path.isfile(join("..", qty, "data", "MD_"+boxName+"_"+str(boxRedshift)+"_"+qty+".fits"))==False :
		hdu2.writeto( join("..", qty, "data", "MD_"+boxName+"_"+str(boxRedshift)+"_"+qty+".fits") )


def convert_pkl_velocity(fileC, fileS, binFile, zList_files, qty):
	"""returns a fits table containing the histograms
	"""
	boxName = fileC.split('/')[5]
	boxZN = float(fileC.split('/')[-1].split('_')[1])
	bins = n.loadtxt(binFile)
	#10**n.arange(0,3.5,0.01)
	dX = ( bins[1:]  - bins[:-1] ) #* n.log(10)
	dlnbin = dX / (( bins[1:]  + bins[:-1] )/2.)

	if boxName=='MD_0.4Gpc' :
		boxLength = 400.
		boxRedshift = 1./boxZN - 1.
		logmp = n.log10(9.63 * 10**7)
		
	if boxName=='MD_1Gpc' :
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
	if os.path.isfile(join("..", qty, "data", "MD_"+boxName+"_"+str(boxRedshift)+"_"+qty+".fits"))==False :
		hdu2.writeto( join("..", qty, "data", "MD_"+boxName+"_"+str(boxRedshift)+"_"+qty+".fits") )


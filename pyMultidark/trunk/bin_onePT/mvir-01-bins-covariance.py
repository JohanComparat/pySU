import glob
from os.path import join
import numpy as n
import astropy.io.fits as fits
import lib_functions_1pt as lib
import os

#Quantity studied
qty = "mvir"

# General information
zList_all =  join(os.environ['PYSU_MD_DIR'], "data", "z-list-all-boxes.txt") 
z0 = n.loadtxt(zList_all,unpack=True)
zList_all2 =  join(os.environ['PYSU_MD_DIR'], "data", "z-list-2LINEAR-COSMO.txt") 
z0short = n.loadtxt(zList_all2,unpack=True,dtype='S')

# redshift lists
dir_boxes =  n.array([os.environ['MD04_DIR'], os.environ['MD10_DIR'], os.environ['MD25_DIR'], os.environ['MD40_DIR'], os.environ['MD25NW_DIR'], os.environ['MD40NW_DIR']])
zList_files = n.array([ join(dir_box,"redshift-list.txt") for dir_box in dir_boxes])

# one point function lists
fileC = n.array(glob.glob( join(os.environ['MULTIDARK_LIGHTCONE_DIR'],"MD_*Gpc*", "properties", qty,"*t_*_Central_JKresampling.pkl")))
fileB = n.array(glob.glob( join( os.environ['MULTIDARK_LIGHTCONE_DIR'],"MD_*Gpc*","properties", qty,"*t_*_"+qty+"_JKresampling.bins")))
fileS = n.array(glob.glob( join( os.environ['MULTIDARK_LIGHTCONE_DIR'],"MD_*Gpc*","properties", qty,"*t_*_Satellite_JKresampling.pkl")))

print "considers ",len(fileC), qty , " function files"

iis = [-1, -2, -3, 3]
for ii in iis:
	lib.plot_CRCoef_mvir(fileC[ii], fileS[ii], fileB[ii],zList_files, z0, z0short, qty)


#rebinned x 2
dataR = n.array([dt[2::2]+dt[1::2] for dt in data])
binsR = bins[1::2]
logmassR = ( binsR[1:]  + binsR[:-1] )/2.
NcountsR = dataR.sum(axis=0) 
okR= ( logmassR> logmp) & (NcountsR>2)

cvR = n.cov(dataR.T[okR])
crR = n.corrcoef(dataR.T[okR])
mmR = logmassR[okR]

mass2XR = interp1d(mmR, n.arange(len(mmR)))

fig = p.figure(0,(6,6))
mat = p.matshow(crR)
p.xticks(n.arange(0,len(mmR),5), mm[n.arange(0,len(mmR),5)],rotation=45)
p.yticks(n.arange(0,len(mmR),5), mm[n.arange(0,len(mmR),5)])
p.axvline(mass2XR(logmp+3), lw=2, color='k')
p.axhline(mass2XR(logmp+3), lw=2, color='k')
cb = p.colorbar(shrink=0.8)
cb.set_label("corrCoef Mvir Hist Counts")
p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
p.ylabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
p.grid()
p.savefig(join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'mvir',"mvir-cr-2.png"))
p.clf()

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

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(log_mvir, log_MF, c=redshift, s=5, marker='o',label="MD "+cos+" data", rasterized=True, vmin=zmin, vmax = zmax)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
p.ylabel(r'log$_{10} (M^2/\rho_m) dn(M)/dM$') 
 # log$_{10}[ n(>M)]')
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)
p.ylim((-4.5,-1))
p.xlim((9.5,16))
p.grid()
p.savefig(join(dir,"mvir-"+figName+cos+"-differential-function-data-xMass.png"))
p.clf()

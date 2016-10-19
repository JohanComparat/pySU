import glob
from os.path import join
import numpy as n
import astropy.io.fits as fits
import lib_functions_1pt as lib
import os
import sys
# data modules
import glob
import sys
import astropy.io.fits as fits
import os
from os.path import join
import cPickle

# numerical modules
import numpy as n
from scipy.interpolate import interp1d
from scipy.misc import derivative
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from scipy.interpolate import interp2d
from scipy.stats import norm


# plotting modules
import matplotlib
#matplotlib.use('pdf')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p

# mass function theory
from hmf import MassFunction
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)


delta_c = 1.686
bh_all = lambda nu, a, b, c : 1+(a**1.5 *nu**2 + a**0.5*b*(a*nu**2)**(1-c) - (a*nu**2)**c/( (a*nu**2)**c+b*(1-c)*(1-c/2)) )/(a**0.5*delta_c)
a=2**(-0.5)
b=0.35
c=0.8
bh = lambda nu : bh_all(nu, a, b, c)
bias = lambda sigma : bh(delta_c/sigma)
"""
xi_mod= lambda R,R0,delta : (R/R0)**(-delta)
xi = lambda R : xi_mod(R, 4, 1.8)

xi(40)

xi(40.)*bh(0.64)*bh(0.749)

xi(400.)*bh(2.85)*bh(3.91)

# nu = delta_c / sigma_M
#ss
"""
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
"""
fileC = n.array(glob.glob( join(os.environ['MULTIDARK_LIGHTCONE_DIR'],"MD_*Gpc*", "properties", qty,"*t_*_Central_JKresampling.pkl")))
fileB = n.array(glob.glob( join( os.environ['MULTIDARK_LIGHTCONE_DIR'],"MD_*Gpc*","properties", qty,"*t_*_"+qty+"_JKresampling.bins")))
fileS = n.array(glob.glob( join( os.environ['MULTIDARK_LIGHTCONE_DIR'],"MD_*Gpc*","properties", qty,"*t_*_Satellite_JKresampling.pkl")))
"""
fileC = n.array(glob.glob( join(os.environ['MULTIDARK_LIGHTCONE_DIR'],"MD_*Gpc*", "properties", qty,"out_*_Central_JKresampling.pkl")))
fileB = n.array(glob.glob( join( os.environ['MULTIDARK_LIGHTCONE_DIR'],"MD_*Gpc*","properties", qty,"out_*_"+qty+"_JKresampling.bins")))
fileS = n.array(glob.glob( join( os.environ['MULTIDARK_LIGHTCONE_DIR'],"MD_*Gpc*","properties", qty,"out_*_Satellite_JKresampling.pkl")))

iis = [8, 13, 29, 31,44]#[-1, -2, -4, -9, -22, 3]
print fileC[iis]


def plot_CRCoef_mvir(fileC, fileS, binFile):
	boxZN = float(os.path.basename(fileC).split('_')[1])
	print boxZN
	hf, boxLength, boxName, boxRedshift, logmp, boxLengthComoving = lib.get_basic_info(fileC, boxZN, delta_wrt='mean')
	bins = n.loadtxt(binFile)
	logmass = ( bins[1:]  + bins[:-1] )/2.
	mass = 10**logmass
	dX = ( 10**bins[1:]  - 10**bins[:-1] )
	dlnbin = (bins[1:]  - bins[:-1])*n.log(10)

	m2sigma = interp1d(hf.M, hf.sigma )
	sigma = m2sigma( mass )
	
	data=cPickle.load(open(fileC,'r'))
	dataS=cPickle.load(open(fileS,'r'))

	cv = n.cov(data.T)
	cr = n.corrcoef(data.T)
	cvS = n.cov(dataS.T)
	crS = n.corrcoef(dataS.T)

	mass2X = interp1d(logmass, n.arange(len(logmass)))
	sigma2X = interp1d(sigma, n.arange(len(logmass)))

	fig = p.figure(0,(6,6))
	mat = p.matshow(cr)
	p.xticks(n.arange(0,len(logmass),5), logmass[n.arange(0,len(logmass),5)],rotation=45)
	p.yticks(n.arange(0,len(logmass),5), logmass[n.arange(0,len(logmass),5)])
	p.axvline(mass2X(logmp+2), lw=2, color='k')
	p.axhline(mass2X(logmp+2), lw=2, color='k')
	#p.axvline(mass2X(logmp+1), lw=2, color='k')
	#p.axhline(mass2X(logmp+1), lw=2, color='k')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("corrCoef")
	p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
	p.ylabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
	p.grid()
	p.savefig(join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'mvir',"covariance","mvir-cr-"+boxName+".png"))
	p.clf()
	
	fig = p.figure(0,(6,6))
	mat = p.matshow(cr)
	p.xticks(n.arange(0,len(sigma),5), n.round(sigma[n.arange(0,len(sigma),5)],3),rotation=45)
	p.yticks(n.arange(0,len(sigma),5), n.round(sigma[n.arange(0,len(sigma),5)],3))
	cb = p.colorbar(shrink=0.8)
	cb.set_label("corrCoef")
	p.xlabel(r'$\sigma$')
	p.ylabel(r'$\sigma$')
	p.grid()
	p.savefig(join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'mvir',"covariance","sigma-cr-"+boxName+".png"))
	p.clf()
	id = int(mass2X(logmp+2.5))
	print id, len(logmass)
	return cr, cv, logmass, sigma, id

crs=n.zeros((len(iis),79,79))
cvs=n.zeros((len(iis),79,79))
lgm=n.zeros((len(iis),79))
sigs=n.zeros((len(iis),79))
ids=n.zeros((len(iis)))
for ii, el in enumerate(iis):
	crs[ii], cvs[ii], lgm[ii], sigs[ii], ids[ii] = plot_CRCoef_mvir(fileC[el], fileS[el], fileB[el])
	
crC=n.zeros((79,79))
crItp=[]
cvC=n.zeros((79,79))

for ii in range(len(iis)):
	crC[ids[ii]:,ids[ii]:] = crs[ii][ids[ii]:,ids[ii]:]
	cvC[ids[ii]:,ids[ii]:] = cvs[ii][ids[ii]:,ids[ii]:]
	xx, yy = n.meshgrid(sigs[ii][ids[ii]:], sigs[ii][ids[ii]:])
	crItp.append(interp2d(xx, yy, crs[ii][ids[ii]:,ids[ii]:], kind='cubic', copy=True, bounds_error=False, fill_value=0.))

logmass = lgm[0]
fig = p.figure(0,(6,6))
mat = p.matshow(crC)
p.xticks(n.arange(0,len(logmass),5), logmass[n.arange(0,len(logmass),5)],rotation=45)
p.yticks(n.arange(0,len(logmass),5), logmass[n.arange(0,len(logmass),5)])
for iid in ids:
	p.axvline(iid-1, lw=1, color='k', ls='dashed')
	p.axhline(iid-1, lw=1, color='k', ls='dashed')

cb = p.colorbar(shrink=0.8)
cb.set_label("corrCoef")
p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
p.ylabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
p.xlim((n.min(ids)-1, 75))
p.ylim((75, n.min(ids)-1))
p.grid()
p.savefig(join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'mvir',"covariance","mvir-cr-all.png"))
p.clf()

#### antidiagonal projection

def getXY(ipos = 40):
	delt = len(logmass) - ipos
	Ntot = len(logmass)
	antiD=[]
	xp=[]
	for ii in n.arange(0, Ntot-ipos, 1):
		antiD.append(crC[ipos+ii][ipos-ii])
		xp.append(logmass[ipos+ii])

	antiD = n.array(antiD)
	xp = n.array(xp)
	p.plot(xp, antiD, label='M='+str(logmass[ipos]))
	return xp, antiD

getXY(30)
getXY(35)
getXY(40)
getXY(45)
getXY(50)
getXY(55)
getXY(60)
getXY(65)
getXY(70)
p.ylabel("corrCoef")
p.xlabel('mass')
p.xlim((logmass[25], 16.))
p.ylim((-0.05, 1.05))
gl = p.legend(loc=0,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig(join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'mvir',"covariance","mvir-antiD.png"))
p.clf()


def getXYh(ipos = 40):
	delt = len(logmass) - ipos
	Ntot = len(logmass)
	antiD=[]
	xp=[]
	for ii in n.arange(0, Ntot-ipos, 1):
		antiD.append(crC[ipos+ii][ipos])
		xp.append(logmass[ipos+ii])

	antiD = n.array(antiD)
	xp = n.array(xp)
	p.plot(xp, antiD, label='M='+str(logmass[ipos]))
	return xp, antiD

getXYh(30)
getXYh(35)
getXYh(40)
getXYh(45)
getXYh(50)
getXYh(55)
getXYh(60)
getXYh(65)
getXYh(70)
p.ylabel("corrCoef")
p.xlabel('mass')
p.xlim((logmass[25], 16.))
p.ylim((-0.05, 1.05))
gl = p.legend(loc=0,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig(join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'mvir',"covariance","mvir-corrHoriZ.png"))
p.clf()

#### SIGMA COVARIANCE

from scipy.interpolate import griddata
# make up some randomly distributed data

xa,ya,za=[],[],[]
for ii in range(len(ids)):
	#x,y position ofinteresting matrix
	xsig_a, ysig_a = n.meshgrid(sigs[ii][ids[ii]:], sigs[ii][ids[ii]:])
	zsig_a = crs[ii][ids[ii]:,ids[ii]:]
	xsig_b = n.ravel(xsig_a)
	ysig_b = n.ravel(ysig_a)
	zsig_b = n.ravel(zsig_a)
	ok = (n.isnan(zsig_b)==False)&(zsig_b!=n.inf)
	xa.append( xsig_b[ok] )
	ya.append( ysig_b[ok] )
	za.append( zsig_b[ok] )

x=n.hstack(xa)
y=n.hstack(ya)
z=n.hstack(za)
	
# define grid.
xi = n.logspace(-0.71, 0.51, 25) #n.arange(0.25, 3.2, 0.02)
yi = n.logspace(-0.71, 0.51, 25) #n.arange(0.25, 3.2, 0.02)
# grid the data.
zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='linear')
# contour the gridded data, plotting dots at the randomly spaced data points.
CS = p.contour(xi,yi,zi,15,linewidths=0.5,colors='k', vmin=0, vmax=1)
CS = p.contourf(xi,yi,zi,15,cmap=p.cm.jet, vmin=0, vmax=1)
cb = p.colorbar(shrink=0.8)
cb.set_label("corrCoef")
p.xlabel(r'$\sigma$')
p.ylabel(r'$\sigma$')
p.plot(x,y,'k,')
p.savefig(join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'mvir',"covariance","griddata.png"))
p.clf()

sss=n.logspace(-1,1,1000)
yyy = bias(sss)**(-0.5)
bmax = n.max(yyy)
bmin = n.min(yyy) 
amin = n.min(zi[n.isnan(zi)==False])
bfun = lambda xx : (bias(xx)**(-0.5)- bmin + amin)/(bmax-bmin+ amin)
xis, yis = n.meshgrid(xi,yi)


zm = bfun(xis) * bfun(yis)
CS = p.contour(xi,yi,zm,15,linewidths=0.5,colors='k', vmin=0, vmax=1)
CS = p.contourf(xi,yi,zm,15,cmap=p.cm.jet, vmin=0, vmax=1)
cb = p.colorbar(shrink=0.8)
cb.set_label(r"$1/(\sqrt{b(\sigma_1)b(\sigma_2)})$")
p.xlabel(r'$\sigma$')
p.ylabel(r'$\sigma$')
p.plot(x,y,'k,')
p.savefig(join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'mvir',"covariance","griddata_model.png"))
p.clf()


CS = p.contour(xi,yi,zi-zm,15,linewidths=0.5,colors='k', vmin=0, vmax=1)
CS = p.contourf(xi,yi,zi-zm,15,cmap=p.cm.jet, vmin=0, vmax=1)
cb = p.colorbar(shrink=0.8)
cb.set_label("residual CC - bias model")
p.xlabel(r'$\sigma$')
p.ylabel(r'$\sigma$')
p.plot(x,y,'k,')
p.savefig(join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'mvir',"covariance","griddata_resid.png"))
p.clf()

tr_val = griddata((x, y), z, (2.6, 2.6), method='linear')-bfun(2.6) * bfun(2.6)
bl_val = griddata((x, y), z, (0.5, 0.5), method='linear')-bfun(0.5) * bfun(0.5)

bnorm = lambda xx : 1-(((bias(xx) - bias(n.min(yi)) )/(bias(n.max(yi))-bias(n.min(yi))) )*(bl_val-tr_val)+tr_val )
 
def res_funs(aa, bb, yVal ): 
	return norm.pdf((bb-aa)/2, loc=0., scale=yVal/20. ) / norm.pdf(0. , loc=0., scale=yVal/20. ) * bnorm(yVal) 

zm2=n.zeros_like(yis)
for jj in range(len(yis)):
	zm2[jj] = res_funs(yis[jj], xis[jj], yi[jj]) 

CS = p.contour(xi,yi,zm2,15,linewidths=0.5,colors='k', vmin=0, vmax=1)
CS = p.contourf(xi,yi,zm2,15,cmap=p.cm.jet, vmin=0, vmax=1)
cb = p.colorbar(shrink=0.8)
cb.set_label("residual")
p.xlabel(r'$\sigma$')
p.ylabel(r'$\sigma$')
p.savefig(join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'mvir',"covariance","griddata_resid_model.png"))
p.clf()


for xid, xx in enumerate(xi):
	ztest = griddata((x, y), z, (xx, yi[yi<xx]), method='linear')
	p.plot(yi[yi<xx],ztest)#,label=str(n.round(xx,2)))

#p.plot(xi, bias(xi)**(-0.5)-(bmax-1),'k--', lw=2, label)
p.plot(xi, bias(xi)**(-0.5)-(bmax-1),'k--', lw=2, label=r"$1+b^{-0.5}_h-max(b^{-0.5}_h)$")
-derivative(bias,xi)
p.ylabel("corrCoef")
p.xlabel('sigma')
p.yscale('log')
p.ylim((0.001, 1.05))
gl = p.legend(loc=0,fontsize=14)
gl.set_frame_on(False)
p.grid()
p.savefig(join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'mvir',"covariance","sigma-corrHoriZ-log.png"))
p.clf()


for xid, xx in enumerate(xi):
	ztest = griddata((x, y), z, (xx, yi[yi<xx]), method='linear')
	p.plot(yi[yi<xx],ztest)#,label=str(n.round(xx,2)))

p.plot(xi, bias(xi)**(-0.5)-(bmax-1),'k--', lw=2, label=r"$1+b^{-0.5}_h-max(b^{-0.5}_h)$")
p.ylabel("corrCoef")
p.xlabel('sigma')
p.ylim((-0.05, 1.05))
gl = p.legend(loc=0,fontsize=14)
gl.set_frame_on(False)
p.grid()
p.savefig(join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'mvir',"covariance","sigma-corrHoriZ.png"))
p.clf()

"""
sss=n.logspace(-1,1,1000)
yyy = bias(sss)**(-0.5)
bmax = n.max(yyy)
p.plot(sss, )
p.savefig(join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'mvir',"covariance","bias-relation.png"))
p.clf()


sigmas = n.arange( 3.5, n.min(sigs), -0.02)
grid_x, grid_y = n.meshgrid(sigmas,sigmas)
#values = 
xx, yy = n.meshgrid(sigs[ii][ids[ii]:], sigs[ii][ids[ii]:])
#crItp.append(interp2d(xx, yy, crs[ii][ids[ii]:,ids[ii]:], 
out = crItp[0]( X,Y)#n.hstack((X)), n.hstack((Y)) )
#mat = out.reshape((len(sigmas),len(sigmas)))

crSS=n.zeros(( len(sigmas), len(sigmas) ))
sigMax = n.array([n.max(sigs[ii][ids[ii]:]) for ii in range(len(iis))])
for ii in range(len(iis)):
	print n.max(sigs[ii][ids[ii]:])
	sel=sigmas[(sigmas<n.max(sigs[ii][ids[ii]:]))]
	idx = n.argmax(sigmas<n.max(sigs[ii][ids[ii]:]))
	crSS[idx:, idx:] = crItp[ii](sel,sel)


fig = p.figure(0,(6,6))
mat = p.matshow(crSS)
p.xticks(n.arange(0,len(sigmas),5), n.round(sigmas[n.arange(0,len(sigmas),5)],3),rotation=45)
p.yticks(n.arange(0,len(sigmas),5), n.round(sigmas[n.arange(0,len(sigmas),5)],3))
cb = p.colorbar(shrink=0.8)
cb.set_label("corrCoef")
p.xlabel(r'$\sigma$')
p.ylabel(r'$\sigma$')
p.xlim((130, 275))
p.ylim((130, 275))
p.grid()
p.savefig(join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'mvir',"covariance","sigma-cr-all.png"))
p.clf()

sys.exit()

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

#rebinning case
dataR = n.array([dt[2::2]+dt[1::2] for dt in data])
binsR = bins[1::2]
logmassR = ( binsR[1:]  + binsR[:-1] )/2.
NcountsR = dataR.sum(axis=0) 
okR= ( logmassR> logmp-0.5) & (NcountsR>2)

cvR = n.cov(dataR.T[okR])
crR = n.corrcoef(dataR.T[okR])
mmR = logmassR[okR]

mass2XR = interp1d(mmR, n.arange(len(mmR)))

fig = p.figure(0,(6,6))
mat = p.matshow(crR)
p.xticks(n.arange(0,len(mmR),5), mmR[n.arange(0,len(mmR),5)],rotation=45)
p.yticks(n.arange(0,len(mmR),5), mmR[n.arange(0,len(mmR),5)])
p.axvline(mass2XR(logmp+3), lw=2, color='k')
p.axhline(mass2XR(logmp+3), lw=2, color='k')
p.axvline(mass2XR(logmp+1), lw=2, color='k')
p.axhline(mass2XR(logmp+1), lw=2, color='k')
cb = p.colorbar(shrink=0.8)
cb.set_label("corrCoef Mvir Hist Counts")
p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
p.ylabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
p.grid()
p.savefig(join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'mvir',"covariance","mvir-cr-2_"+figName+".png"))
p.clf()



Ncounts = data.sum(axis=0) 
Nall = Ncounts / volume
ok= ( sigma> logmp-0.5) & (Ncounts>2)

index=n.arange(int(data.shape[0]))
n.random.shuffle( index )
Ntotal = int(data.shape[0])

dataS = n.array([n.sum(data[id:id+Ntotal/resamp:1], axis=0) for id in n.arange(0,Ntotal,Ntotal/resamp)])

cv = n.cov(data.T[ok])
cr = n.corrcoef(data.T[ok])
mm = logmass[ok]
sigma = sig[ok]
nu = nus[ok]

cvS = n.cov(dataS.T[ok])
crS = n.corrcoef(dataS.T[ok])

mass2X = interp1d(mm, n.arange(len(mm)))

fig = p.figure(0,(6,6))
mat = p.matshow(cr)
p.xticks(n.arange(0,len(mm),5), mm[n.arange(0,len(mm),5)],rotation=45)
p.yticks(n.arange(0,len(mm),5), mm[n.arange(0,len(mm),5)])
p.axvline(mass2X(logmp+3), lw=2, color='k')
p.axhline(mass2X(logmp+3), lw=2, color='k')
p.axvline(mass2X(logmp+1), lw=2, color='k')
p.axhline(mass2X(logmp+1), lw=2, color='k')
cb = p.colorbar(shrink=0.8)
cb.set_label("corrCoef Mvir Counts "+figName)
p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
p.ylabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
p.grid()
p.savefig(join(os.environ['MULTIDARK_LIGHTCONE_DIR'], 'mvir',"covariance","mvir-cr-0_"+figName+".png"))
p.clf()

fig = p.figure(0,(6,6))
mat = p.matshow(cv)
p.xticks(n.arange(0,len(nu),5), n.round(nu[n.arange(0,len(nu),5)],3),rotation=45)
p.yticks(n.arange(0,len(nu),5), n.round(nu[n.arange(0,len(nu),5)],3))
#p.axvline(mass2X(logmp+3), lw=2, color='k')
#p.axhline(mass2X(logmp+3), lw=2, color='k')
#p.axvline(mass2X(logmp+1), lw=2, color='k')
#p.axhline(mass2X(logmp+1), lw=2, color='k')
cb = p.colorbar(shrink=0.8)
cb.set_label("corrCoef Mvir Counts "+figName)
p.xlabel(r'$\nu$')
p.ylabel(r'$\nu$')
p.grid()

"""

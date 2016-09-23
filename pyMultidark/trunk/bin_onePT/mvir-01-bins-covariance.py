import glob
from os.path import join
import numpy as n
import astropy.io.fits as fits
import lib_functions_1pt as lib
import os
import sys

delta_c = 1.686
bh_all = lambda nu, a, b, c : 1+(a**1.5 *nu**2 + a**0.5*b*(a*nu**2)**(1-c) - (a*nu**2)**c/( (a*nu**2)**c+b*(1-c)*(1-c/2)) )/(a**0.5*delta_c)
a=2**(-0.5)
b=0.35
c=0.8
bh = lambda nu : bh_all(nu, a, b, c)

xi_mod= lambda R,R0,delta : (R/R0)**(-delta)
xi = lambda R : xi_mod(R, 4, 1.8)

# nu = delta_c / sigma_M
#ss

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

iis = [-1, -2, -4, -9, -22, 3]
for ii in iis:
	mm, sigma, nu = lib.plot_CRCoef_mvir(fileC[ii], fileS[ii], fileB[ii],zList_files, z0, z0short, qty,rebin=False)
	print mm, nu

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

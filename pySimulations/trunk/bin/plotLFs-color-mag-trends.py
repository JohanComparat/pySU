#! /usr/bin/env python

"""
This script produces quality plots to check that the LFs are fine compared to simumlations.
"""

import sys
import os
from os.path import join
data_dir = os.environ['DATA_DIR']
import glob

from lib_plot import *
#from lineListAir import *
SNlim = 5


plotDir="/home/comparat/database/Simulations/galform-lightcone/products/emissionLineLuminosityFunctions/plots/"

dir="/home/comparat/database/Simulations/galform-lightcone/products/emissionLineLuminosityFunctions/O2_3728/"

"VVDSDEEP-MagLimI-22.5"
"DEEP2-MagLimR-"



lf_measurement_files=n.array(glob.glob(dir+"*MagLimI-*.txt"))
lf_measurement_files.sort()


lf_ref=lf_measurement_files[-1]
dataRef = n.loadtxt( lf_ref, unpack=True)
phiRatio = n.empty([ len(lf_measurement_files[:-1]), len(dataRef[0]) ])
label = n.array(["I<22.0", "I<22.5", "I<23.0","I<23.5"])
for jj in range(len(lf_measurement_files[:-1])):
	lf_obs=lf_measurement_files[jj]
	data= n.loadtxt( lf_obs, unpack=True)
	phiRatio[jj] = data[3] / dataRef[3]	
	print lf_obs
	

p.figure(0,(6,6))
for jj in range(len(label)):
	p.plot(dataRef[2],phiRatio[jj],label=label[jj])

p.xlabel(r'$log_{10}(L[O_{II}])$ [erg s$^{-1}$]')
p.ylabel(r'$\Phi/\Phi_{ref}$')
p.xscale('log')
p.xlim((1e40,1e43))
p.ylim((-0.05,1.05))
p.grid()
p.legend(loc=2)
p.savefig(join(plotDir,"O2_3728_VVDSDEEP_trendMag-z0.7.jpg"))
p.clf()



#####################################3
#####################################3
# R-Z
#####################################3
#####################################3

lf_measurement_files_ref=n.array(glob.glob(dir+"*VVDSrz_gt_0.0-z0.7*.txt"))
lf_measurement_files=n.array(glob.glob(dir+"*VVDSrz_?t_*z0.7*.txt"))
lf_measurement_files.sort()

dataRef = n.loadtxt( lf_measurement_files_ref[0], unpack=True)

label = n.array(["r-z>0", "r-z>0.5", "r-z>1", "r-z>1.5", "r-z<1", "r-z<1.5", "r-z<2"])
phiRatio = n.empty([ 7, len(dataRef[0]) ])
for ii, el in enumerate(lf_measurement_files):
	data= n.loadtxt( el, unpack=True)
	phiRatio[ii] = data[3] / dataRef[3]	

p.figure(0,(6,6))
for jj in range(len(label)):
	p.plot(dataRef[2],phiRatio[jj],label=label[jj])

p.xlabel(r'$log_{10}(L[O_{II}])$ [erg s$^{-1}$]')
p.ylabel(r'$\Phi/\Phi_{ref}$')
p.xscale('log')
p.xlim((7e40,5e43))
p.ylim((-0.05,1.05))
p.grid()
p.legend(loc=4)
p.savefig(join(plotDir,"trends_O2_3728_I22.5_RZ-z0.75.pdf"))
p.clf()


lf_measurement_files_ref=n.array(glob.glob(dir+"*VVDSrz_gt_0.0-z0.9*.txt"))
lf_measurement_files=n.array(glob.glob(dir+"*VVDSrz_?t_*z0.9*.txt"))
lf_measurement_files.sort()

dataRef = n.loadtxt( lf_measurement_files_ref[0], unpack=True)

phiRatio = n.empty([ 7, len(dataRef[0]) ])
for ii, el in enumerate(lf_measurement_files):
	data= n.loadtxt( el, unpack=True)
	phiRatio[ii] = data[3] / dataRef[3]	

p.figure(0,(6,6))
for jj in range(len(label)):
	p.plot(dataRef[2],phiRatio[jj],label=label[jj])

p.xlabel(r'$log_{10}(L[O_{II}])$ [erg s$^{-1}$]')
p.ylabel(r'$\Phi/\Phi_{ref}$')
p.xscale('log')
p.xlim((7e40,5e43))
p.ylim((-0.05,1.05))
p.grid()
p.legend(loc=4)
p.savefig(join(plotDir,"trends_O2_3728_I22.5_RZ-z0.9.pdf"))
p.clf()


lf_measurement_files_ref=n.array(glob.glob(dir+"*VVDSrz_gt_0.0-z1.*.txt"))
lf_measurement_files=n.array(glob.glob(dir+"*VVDSrz_?t_*z1.*.txt"))
lf_measurement_files.sort()

dataRef = n.loadtxt( lf_measurement_files_ref[0], unpack=True)

phiRatio = n.empty([ 7, len(dataRef[0]) ])
for ii, el in enumerate(lf_measurement_files):
	data= n.loadtxt( el, unpack=True)
	phiRatio[ii] = data[3] / dataRef[3]	

p.figure(0,(6,6))
for jj in range(len(label)):
	p.plot(dataRef[2],phiRatio[jj],label=label[jj])

p.xlabel(r'$log_{10}(L[O_{II}])$ [erg s$^{-1}$]')
p.ylabel(r'$\Phi/\Phi_{ref}$')
p.xscale('log')
p.xlim((7e40,5e43))
p.ylim((-0.05,1.05))
p.grid()
p.legend(loc=4)
p.savefig(join(plotDir,"trends_O2_3728_I22.5_RZ-z1.2.pdf"))
p.clf()



#####################################3
#####################################3
# G-R
#####################################3
#####################################3


lf_measurement_files_ref=n.array(glob.glob(dir+"*VVDSgr_gt_0.0-z0.7*.txt"))
lf_measurement_files=n.array(glob.glob(dir+"*VVDSgr_?t_*z0.7*.txt"))
lf_measurement_files.sort()

dataRef = n.loadtxt( lf_measurement_files_ref[0], unpack=True)

label = n.array(["g-r>0", "g-r>0.5", "g-r>1", "g-r>1.5", "g-r<1", "g-r<1.5", "g-r<2"])
phiRatio = n.empty([ 7, len(dataRef[0]) ])
for ii, el in enumerate(lf_measurement_files):
	data= n.loadtxt( el, unpack=True)
	phiRatio[ii] = data[3] / dataRef[3]	

p.figure(0,(6,6))
for jj in range(len(label)):
	p.plot(dataRef[2],phiRatio[jj],label=label[jj])

p.xlabel(r'$log_{10}(L[O_{II}])$ [erg s$^{-1}$]')
p.ylabel(r'$\Phi/\Phi_{ref}$')
p.xscale('log')
p.xlim((7e40,5e43))
p.ylim((-0.05,1.05))
p.grid()
p.legend(loc=4)
p.savefig(join(plotDir,"trends_O2_3728_I22.5_GR-z0.75.pdf"))
p.clf()


lf_measurement_files_ref=n.array(glob.glob(dir+"*VVDSgr_gt_0.0-z0.9*.txt"))
lf_measurement_files=n.array(glob.glob(dir+"*VVDSgr_?t_*z0.9*.txt"))
lf_measurement_files.sort()

dataRef = n.loadtxt( lf_measurement_files_ref[0], unpack=True)

label = n.array(["g-r>0", "g-r>0.5", "g-r>1", "g-r>1.5", "g-r<1", "g-r<1.5", "g-r<2"])
phiRatio = n.empty([ 7, len(dataRef[0]) ])
for ii, el in enumerate(lf_measurement_files):
	data= n.loadtxt( el, unpack=True)
	phiRatio[ii] = data[3] / dataRef[3]	

p.figure(0,(6,6))
for jj in range(len(label)):
	p.plot(dataRef[2],phiRatio[jj],label=label[jj])

p.xlabel(r'$log_{10}(L[O_{II}])$ [erg s$^{-1}$]')
p.ylabel(r'$\Phi/\Phi_{ref}$')
p.xscale('log')
p.xlim((7e40,5e43))
p.ylim((-0.05,1.05))
p.grid()
p.legend(loc=4)
p.savefig(join(plotDir,"trends_O2_3728_I22.5_GR-z0.9.pdf"))
p.clf()


lf_measurement_files_ref=n.array(glob.glob(dir+"*VVDSgr_gt_0.0-z1.*.txt"))
lf_measurement_files=n.array(glob.glob(dir+"*VVDSgr_?t_*z1.*.txt"))
lf_measurement_files.sort()

dataRef = n.loadtxt( lf_measurement_files_ref[0], unpack=True)

label = n.array(["g-r>0", "g-r>0.5", "g-r>1", "g-r>1.5", "g-r<1", "g-r<1.5", "g-r<2"])
phiRatio = n.empty([ 7, len(dataRef[0]) ])
for ii, el in enumerate(lf_measurement_files):
	data= n.loadtxt( el, unpack=True)
	phiRatio[ii] = data[3] / dataRef[3]	

p.figure(0,(6,6))
for jj in range(len(label)):
	p.plot(dataRef[2],phiRatio[jj],label=label[jj])

p.xlabel(r'$log_{10}(L[O_{II}])$ [erg s$^{-1}$]')
p.ylabel(r'$\Phi/\Phi_{ref}$')
p.xscale('log')
p.xlim((7e40,5e43))
p.ylim((-0.05,1.05))
p.grid()
p.legend(loc=4)
p.savefig(join(plotDir,"trends_O2_3728_I22.5_GR-z1.2.pdf"))
p.clf()


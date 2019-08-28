#! /usr/bin/env python

"""
This script produces the stacks for emission line luminosity limited samples.

nohup python2 stack_spectra_ELG_LineLF.py > stack_spectra_ELG_LineLF.log &

"""
import sys
import os 
from os.path import join
import glob
import numpy as n
import SpectraStackingEBOSS as sse
from scipy.interpolate import interp1d

spec_dir = join(os.environ['HOME'],"SDSS/stacks/AW/")
"""
spec_dir2 = join(os.environ['HOME'],"SDSS/stacks/v0_LineLF/")

PLATE, MJD, FIBREID, z, z_err, log_l_O2, log_l_O3, log_l_Hb, WEIGHT_SYSTOT, WEIGHT_CP, WEIGHT_NOZ = n.loadtxt(os.path.join(spec_dir2, "ELGv5_11_0rrv2_all_lum_weight.txt"), unpack=True)

weight=WEIGHT_SYSTOT*WEIGHT_CP* WEIGHT_NOZ

data_to_stack = n.transpose([PLATE, MJD, FIBREID, z, weight])

both = (log_l_O2>0)&(log_l_O3>0)
deltaO2O3 = log_l_O2 - log_l_O3

name = 'elg_w_O2O3_zmin_0.5_zmax_1.05.ascii'
print(name)
z_sel = (z>=0.5)&(z<1.05)&(both)&(weight>0)
n.savetxt( os.path.join(spec_dir, name),data_to_stack[z_sel])

# per redshift bins
outN, outBins = n.histogram(z[(z>0.5) & (z<1.4) & (both)], n.arange(z.min(),z.max(), 0.001))

outNC = n.cumsum(outN)

itp = interp1d(outNC, outBins[1:])

zmins = itp(n.arange(200,len(z[both]),10000))

for zmin, zmax in zip(zmins[:-1], zmins[1:]):
	z_sel = (z>=zmin)&(z<zmax)&(both)&(weight>0)
	print(zmin, zmax, len(z[z_sel]), n.min(z[z_sel]), n.max(z[z_sel]))
	name = 'elg_w_O2O3_zmin_'+str(n.round(zmin,4))+'_zmax_'+str(n.round(zmax,4))+'.ascii'
	print(name)
	n.savetxt( os.path.join(spec_dir, name),data_to_stack[z_sel])

# per deltaO2O3 bins
outN, outBins = n.histogram(deltaO2O3[(z>0.5) & (z<1.4) & (both)], n.arange(n.min(deltaO2O3[(z>0.5) & (z<1.4) & (both)]),n.max(deltaO2O3[(z>0.5) & (z<1.4) & (both)]), 0.001))

outNC = n.cumsum(outN)

itp = interp1d(outNC, outBins[1:])

zmins = itp(n.arange(200,len(z[both]),10000))

for zmin, zmax in zip(zmins[:-1], zmins[1:]):
	z_sel = (deltaO2O3>=zmin)&(deltaO2O3<zmax)&(both)&(weight>0)
	print(zmin, zmax, len(z[z_sel]), n.min(deltaO2O3[z_sel]), n.max(deltaO2O3[z_sel]))
	name = 'elg_w_O2O3_O2O3min_'+str(n.round(zmin,4))+'_O2O3max_'+str(n.round(zmax,4))+'.ascii'
	print(name)
	n.savetxt( os.path.join(spec_dir, name),data_to_stack[z_sel])


"""

file_list = n.array(glob.glob(os.path.join(spec_dir, '*.ascii')))

def stack_it(specList ):
	outfile = join(spec_dir, os.path.basename(specList)[:-6]+".stack")
	stack=sse.SpectraStackingEBOSS(specList, outfile, PBKT_input=True   )
	stack.createStackMatrix()
	print(outfile)
	stack.stackSpectra()

for file_input in file_list:
	stack_it(file_input)


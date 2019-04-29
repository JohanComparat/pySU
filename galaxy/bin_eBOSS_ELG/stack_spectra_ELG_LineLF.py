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

spec_dir = join(os.environ['HOME'],"SDSS/stacks/LineLF/")
"""
PLATE, MJD, FIBREID, z, z_err, log_l_O2, log_l_O3, log_l_Hb, WEIGHT_SYSTOT, WEIGHT_CP, WEIGHT_NOZ = n.loadtxt(os.path.join(spec_dir, "ELGv5_11_0rrv2_all_lum_weight.txt"), unpack=True)

weight=WEIGHT_SYSTOT*WEIGHT_CP* WEIGHT_NOZ

N_per_stack=1000

zmins = n.arange(0.6,1.2,0.2)
zmaxs = zmins+0.2

for zmin, zmax in zip(zmins, zmaxs):
  z_sel = (z>=zmin)&(z<zmax)&(log_l_O2>0)&(weight>0)
  idx_L = n.argsort(log_l_O2[z_sel])
  #log_l_O2[z_sel][idx_L]
  data_to_stack = n.transpose([PLATE, MJD, FIBREID, z, weight])[z_sel][idx_L]
  # write in files of 1000
  ids_bds = n.arange(0, len(data_to_stack), N_per_stack)
  for ids_bd in ids_bds:
    name = 'elg_'+str(n.min(log_l_O2[z_sel][idx_L][::-1][ids_bd:ids_bd+N_per_stack])).zfill(3)+'_O2_'+str(n.max(log_l_O2[z_sel][idx_L][::-1][ids_bd:ids_bd+N_per_stack])).zfill(3)+"-"+str(zmin)+'_z_'+str(zmax)+'_N_'+str(len(log_l_O2[z_sel][idx_L][ids_bd:ids_bd+N_per_stack]))+'.ascii'
    print(name)
    n.savetxt( os.path.join(spec_dir, name),data_to_stack[::-1][ids_bd:ids_bd+N_per_stack])


zmins = n.arange(0.6,0.8,0.2)
zmaxs = zmins+0.2

for zmin, zmax in zip(zmins, zmaxs):
  z_sel = (z>=zmin)&(z<zmax)&(log_l_O3>0)&(weight>0)
  idx_L = n.argsort(log_l_O3[z_sel])
  #log_l_O3[z_sel][idx_L]
  data_to_stack = n.transpose([PLATE, MJD, FIBREID, z, weight])[z_sel][idx_L]
  # write in files of 1000
  ids_bds = n.arange(0, len(data_to_stack), N_per_stack)
  for ids_bd in ids_bds:
    name = 'elg_'+str(n.min(log_l_O3[z_sel][idx_L][::-1][ids_bd:ids_bd+N_per_stack])).zfill(3)+'_O3_'+str(n.max(log_l_O3[z_sel][idx_L][::-1][ids_bd:ids_bd+N_per_stack])).zfill(3)+"-"+str(zmin)+'_z_'+str(zmax)+'_N_'+str(len(log_l_O3[z_sel][idx_L][ids_bd:ids_bd+N_per_stack]))+'.ascii'
    print(name)
    n.savetxt( os.path.join(spec_dir, name),data_to_stack[::-1][ids_bd:ids_bd+N_per_stack])


zmins = n.arange(0.6,0.8,0.2)
zmaxs = zmins+0.2

for zmin, zmax in zip(zmins, zmaxs):
  z_sel = (z>=zmin)&(z<zmax)&(log_l_Hb>0)&(weight>0)
  idx_L = n.argsort(log_l_Hb[z_sel])
  #log_l_Hb[z_sel][idx_L]
  data_to_stack = n.transpose([PLATE, MJD, FIBREID, z, weight])[z_sel][idx_L]
  # write in files of 1000
  ids_bds = n.arange(0, len(data_to_stack), N_per_stack)
  for ids_bd in ids_bds:
    name = 'elg_'+str(n.min(log_l_Hb[z_sel][idx_L][::-1][ids_bd:ids_bd+N_per_stack])).zfill(3)+'_Hb_'+str(n.max(log_l_Hb[z_sel][idx_L][::-1][ids_bd:ids_bd+N_per_stack])).zfill(3)+"-"+str(zmin)+'_z_'+str(zmax)+'_N_'+str(len(log_l_Hb[z_sel][idx_L][ids_bd:ids_bd+N_per_stack]))+'.ascii'
    print(name)
    n.savetxt( os.path.join(spec_dir, name),data_to_stack[::-1][ids_bd:ids_bd+N_per_stack])
"""

file_list = n.array(glob.glob(os.path.join(spec_dir, '*.ascii')))

def stack_it(specList ):
	outfile = join(spec_dir, os.path.basename(specList)[:-6]+".stack")
	stack=sse.SpectraStackingEBOSS(specList, outfile, PBKT_input=True   )
	stack.createStackMatrix_Weighted()
	print(outfile)
	stack.stackSpectra()

for file_input in file_list:
	stack_it(file_input)


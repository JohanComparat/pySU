#! /usr/bin/env python

"""
This script produces the stacks for samples defined by a list of spectra identifiers.

nohup python3 stack_spectra_ELG.py > stack_spectra_ELG.log &

"""
import sys
import os 
from os.path import join
import glob
import numpy as n
import SpectraStackingEBOSS as sse

# T2 templates
~/SDSS/stacks/X_AGN$ ls -lh ROSAT_AGNT2*stack
ROSAT_AGNT2_zmin_00_zmax_02..stack for 6100 AA < lambda < 8100 AA
ROSAT_AGNT2_zmin_00_zmax_05..stack for 3650 AA < lambda < 6100 AA
ROSAT_AGNT2_zmin_03_zmax_08..stack for 3000 AA < lambda < 3650 AA
ROSAT_AGNT2_zmin_05_zmax_10..stack for 2300 AA < lambda < 3000 AA

adjust heights of each template so they match in the boundary region

Write a template in the 4MOST format.

# T1 template
~/SDSS/stacks/X_AGN
ROSAT_AGNT1_zmin_*..stack

stack_dir = join(os.environ['HOME'],"SDSS/stacks")

def stack_it( specList ):
	outfile = join(stack_dir, os.path.basename(specList)+".stack")
	print(stack_dir, outfile)
	stack=sse.SpectraStackingEBOSS(specList, outfile )
	stack.createStackMatrix()
	stack.stackSpectra()

list_2_stack = n.array(glob.glob(join(stack_dir, "*.ascii")))
for el in list_2_stack:
	stack_it(el)

for el in list_2_stack[::-1]:
	stack_it(el)

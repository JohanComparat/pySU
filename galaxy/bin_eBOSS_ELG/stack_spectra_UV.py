#! /usr/bin/env python

"""
This script produces the stacks for emission line luminosity limited samples.
"""
import sys
import os 
from os.path import join
import glob
import numpy as n
import SpectraStackingEBOSS as sse
import time
t0=time.time()
# create all stacks
dataList = n.array(glob.glob(join(os.environ['HOME'],"SDSS/lss/catalogs/3/stacks_v1", "eboss-elg_*.asc")))
dataList.sort()

for specList in dataList:
	print('considers', specList, time.time()-t0)
	outfile = join(os.environ['HOME'], "stacks", os.path.basename(specList)[:-4]+".stack")
	#outfile = join(os.environ['HOME'], "SDSS", "stacks", os.path.basename(specList)[:-4]+".stack")
	#if os.path.isfile(outfile)==False:
	print('starts working on', outfile, time.time()-t0)
	stack=sse.SpectraStackingEBOSS(specList, outfile)
	if os.path.isfile(outfile+'.specMatrix.dat')==False:
		print('creates matrix', time.time()-t0)
		stack.createStackMatrix()
	print('stacks', time.time()-t0)
	stack.stackSpectra()

	outfile = join(os.environ['HOME'], "stacks", os.path.basename(specList)[:-4]+".UVstack")
	#outfile = join(os.environ['HOME'],"SDSS", "stacks", os.path.basename(specList)[:-4]+".UVstack")
	#if os.path.isfile(outfile)==False:
	print('starts working on', outfile, time.time()-t0)
	stack=sse.SpectraStackingEBOSS(specList, outfile)
	if os.path.isfile(outfile+'.specMatrix.dat')==False:
		print('creates matrix', time.time()-t0)
		stack.createStackMatrix_UVnormed()
	print('stacks', time.time()-t0)
	stack.stackSpectra()
		

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

# create all stacks
dataList = n.array(glob.glob(join(os.environ['HOME'],"SDSS/lss/catalogs/3/stacks_v1", "eboss-elg_*.asc")))
dataList.sort()

for specList in dataList:
	print( specList )
	outfile = join(os.environ['HOME'], "stacks", os.path.basename(specList)[:-4]+".stack")
	#outfile = join(os.environ['HOME'], "SDSS", "stacks", os.path.basename(specList)[:-4]+".stack")
	if os.path.isfile(outfile)==False:
		stack=sse.SpectraStackingEBOSS(specList, outfile)
		print(outfile)
		if os.path.isfile(outfile+'.specMatrix.dat')==False:
			stack.createStackMatrix()
		stack.stackSpectra()

	outfile = join(os.environ['HOME'], "stacks", os.path.basename(specList)[:-4]+".UVstack")
	#outfile = join(os.environ['HOME'],"SDSS", "stacks", os.path.basename(specList)[:-4]+".UVstack")
	if os.path.isfile(outfile)==False:
		stack=sse.SpectraStackingEBOSS(specList, outfile)
		print(outfile)
		if os.path.isfile(outfile+'.specMatrix.dat')==False:
			stack.createStackMatrix_UVnormed()
		stack.stackSpectra()
		

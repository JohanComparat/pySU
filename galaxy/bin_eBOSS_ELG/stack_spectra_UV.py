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
dataList = n.array(glob.glob(join(os.environ['HOME'],"SDSS/stacks", "eboss-elg_*.asc")))
dataList.sort()

for specList in dataList:
	#specList = dataList[-1]
	print('considers', specList, time.time()-t0)
	outfile = join(os.environ['HOME'], "SDSS/stacks", os.path.basename(specList)[:-4]+".stack")
	#outfile = join(os.environ['HOME'], "SDSS", "stacks", os.path.basename(specList)[:-4]+".stack")
	#if os.path.isfile(outfile)==False:
	print('starts working on', outfile, time.time()-t0)
	stack=sse.SpectraStackingEBOSS(specList, outfile, dLambda = 0.0002, dV=-9999.99, l_start=3.35, l_end=3.579)
	if os.path.isfile(outfile+'.specMatrix.dat')==False:
		print('creates matrix', time.time()-t0)
		stack.createStackMatrix()
	print('stacks', time.time()-t0)
	stack.stackSpectra()

hdu = fits.open(outfile)
dd=hdu[1].data
wl=dd['wavelength'         ]
s1 = (dd['NspectraPerPixel'   ]>0.5*n.max(dd['NspectraPerPixel'   ]))
x, y= dd['wavelength'         ][s1], dd['medianStack'        ][s1]
yerr = dd['jackknifStackErrors'][s1]
pfit = stack.fit_UV_continuum(x,y,yerr,degree=3)
Fcont = n.polyval(pfit, dd['wavelength'         ])
spec_UV.append(dd['medianStack'        ]/Fcont)
spec_UV_Err.append(dd['jackknifStackErrors'        ]/Fcont)



outfile = join(os.environ['HOME'], "SDSS/stacks", os.path.basename(specList)[:-4]+".UVstack")


if os.path.isfile(self.out_file):
	os.remove(self.out_file)
print( "stack written to", self.out_file )
thdulist.writeto(self.out_file)


	outfile = join(os.environ['HOME'], "SDSS/stacks", os.path.basename(specList)[:-4]+".UVstack")
	#outfile = join(os.environ['HOME'],"SDSS", "stacks", os.path.basename(specList)[:-4]+".UVstack")
	#if os.path.isfile(outfile)==False:
	print('starts working on', outfile, time.time()-t0)
	stack=sse.SpectraStackingEBOSS(specList, outfile, dLambda = 0.0002, dV=-9999.99, l_start=3.34, l_end=3.556)
	if os.path.isfile(outfile+'.specMatrix.dat')==False:
		print('creates matrix', time.time()-t0)
		stack.createStackMatrix_UVnormed()
	print('stacks', time.time()-t0)
	stack.stackSpectra()

#specMatrix = n.loadtxt(outfile+'.specMatrix.dat')
#specMatrixErr = n.loadtxt(outfile+'.specMatrixErr.dat')
#specMatrixWeight = n.loadtxt(outfile+'.specMatrixWeight.dat')

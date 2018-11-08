#! /usr/bin/env python

"""
This script produces the stacks for emission line luminosity limited samples.
"""
import astropy.io.fits as fits
import sys
import os 
from os.path import join
import glob
import numpy as n
import SpectraStackingEBOSS as sse
import time
t0=time.time()
# create all stacks
#stack_dir = join( os.environ['HOME'], "SDSS/stacks/v2" )
stack_dir = join( os.environ['HOME'], "SDSS/stacks" )
dataList = n.array(glob.glob(join(stack_dir, "eboss-elg_*.asc")))
dataList.sort()

for specList in dataList:
	print('considers', specList, time.time()-t0)
	out_file = join(stack_dir, os.path.basename(specList)[:-4]+".stack")
	out_file_UV = join(stack_dir, os.path.basename(specList)[:-4]+".UVstack")
	print('starts working on', out_file, time.time()-t0)
	stack=sse.SpectraStackingEBOSS(specList, out_file, dLambda = 0.0001, dV=-9999.99, l_start=3.35, l_end=3.579)
	if os.path.isfile(out_file+'.specMatrix.dat')==False:
		print('creates matrix', time.time()-t0)
		stack.createStackMatrix()
	if os.path.isfile(out_file)==False:
		print('stacks', time.time()-t0)
		stack.stackSpectra()
	hdu = fits.open(out_file)
	dd=hdu[1].data
	wl=dd['wavelength'         ]
	s1 = (dd['NspectraPerPixel'   ]>0.5*n.max(dd['NspectraPerPixel'   ]))
	x, y= dd['wavelength'         ][s1], dd['medianStack'        ][s1]
	yerr = dd['jackknifStackErrors'][s1]
	pfit = stack.fit_UV_continuum(x,y,yerr,degree=5)
	Fcont = n.polyval(pfit, dd['wavelength'         ])
	spec_UV = dd['medianStack'        ]/Fcont
	spec_UV_Err = dd['jackknifStackErrors'        ]/Fcont
	c1 = fits.Column(name="medianStack_UVnormed",format="D", unit="erg/s/cm2/Angstrom", array= spec_UV)
	c2 = fits.Column(name="jackknifStackErrors_UVnormed",format="D", unit="erg/s/cm2/Angstrom", array= spec_UV_Err)
	cols_2 = fits.ColDefs([c1, c2])
	tbhdu = fits.BinTableHDU.from_columns(hdu[1].columns + cols_2)
	prihdr = fits.Header()
	prihdr['author'] = "JC"
	prihdr['survey'] = hdu[0].header['SURVEY']
	prihdr['in_file'] = os.path.basename(specList)[:-4]
	prihdr['Nspec'] = hdu[0].header['Nspec']
	prihdu = fits.PrimaryHDU(header=prihdr)
	thdulist = fits.HDUList([prihdu, tbhdu])
	if os.path.isfile(out_file_UV):
		os.remove(out_file_UV)
	print( "stack written to", out_file_UV )
	thdulist.writeto(out_file_UV)


#specMatrix = n.loadtxt(out_file+'.specMatrix.dat')
#specMatrixErr = n.loadtxt(out_file+'.specMatrixErr.dat')
#specMatrixWeight = n.loadtxt(out_file+'.specMatrixWeight.dat')

#out_file = join(os.environ['HOME'], "SDSS/stacks", os.path.basename(specList)[:-4]+".UVstack")
##out_file = join(os.environ['HOME'],"SDSS", "stacks", os.path.basename(specList)[:-4]+".UVstack")
##if os.path.isfile(out_file)==False:
#print('starts working on', out_file, time.time()-t0)
#stack=sse.SpectraStackingEBOSS(specList, out_file, dLambda = 0.0002, dV=-9999.99, l_start=3.34, l_end=3.556)
#if os.path.isfile(out_file+'.specMatrix.dat')==False:
	#print('creates matrix', time.time()-t0)
	#stack.createStackMatrix_UVnormed()
#print('stacks', time.time()-t0)
#stack.stackSpectra()

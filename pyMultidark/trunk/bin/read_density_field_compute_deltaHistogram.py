#  cd pySU/pyMultidark/trunk/bin/fortranfile-0.2.1/

import numpy as n
import os
from os.path import join
from astropy.io import fits
import time
import fortranfile
DFdir = join("/data2", "users", "gustavo", "BigMD", "1Gpc_3840_Planck1_New", "DENSFIELDS")
mockDir = "/data1/DATA/eBOSS/Multidark-box-mocks/parts/"

inFiles = n.array(["dmdens_cic_104.dat", "dmdens_cic_101.dat", "dmdens_cic_097.dat", "dmdens_cic_087.dat"])

for infi in inFiles:
	bins = n.hstack((0,n.logspace(-3, 4, 10000)))
	DFfile = join(DFdir,infi)
	f = fortranfile.FortranFile(DFfile)
	gridx, gridy, gridz = f.readInts()
	res = n.empty((gridx, len(bins)-1))
	for kk in range(gridx):
		DF = f.readReals()
		res[kk],bi = n.histogram(DF,bins=bins)

	result = n.sum(res,axis=0)
	path_to_outputCat =  join(mockDir,infi[:-4] + "_DFhist.dat")
	f=open(path_to_outputCat, 'w')
	cPickle.dump([bins,result],f)
	f.close()

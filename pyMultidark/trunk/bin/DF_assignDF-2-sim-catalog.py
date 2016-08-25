#  cd pySU/pyMultidark/trunk/bin/fortranfile-0.2.1/

import numpy as n
import os
from os.path import join
from astropy.io import fits
import time
import fortranfile
import glob

# loads the density field
DFdir = join("/data2", "users", "gustavo", "BigMD", "1Gpc_3840_Planck1_New", "DENSFIELDS")
DFfile = join(DFdir,"dmdens_cic_087.dat")

# loads the snapshot files 
snList= n.array(glob.glob( "/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.403200_PM_Nb_?.fits"))#
def writeDFMock(dataCat, DFfile, Lbox = 1000.):
	print dataCat, DFfile
	md = fits.open(dataCat)[1].data
	path_to_outputCat = dataCat[:-4] + "DF.fits.gz"
	# opens the DF file
	f = fortranfile.FortranFile(DFfile)
	gridx, gridy, gridz = f.readInts()
	dx = Lbox/gridx
	# convert QSO positions into indexes
	i = ( ( md['x'] / gridx ) // 1 ).astype( 'int' )
	j = ( ( md['y'] / gridx ) // 1 ).astype( 'int' )
	k= ( ( md['z'] / gridx ) // 1 ).astype( 'int' )
	#init the output array :
	delta = n.ones_like(i)*-1.
	#delta1 = n.empty_like(x)
	#delta2 = n.empty_like(x)
	# now loops over k (every line of the file) and assigns delta values.
	for kk in range(gridx):
		print kk
		sel = (k==kk)
		N = i[sel] + gridx * j[sel] 
		DF = f.readReals()
		delta[sel] = DF[N]
		
	f.close()

	c0 = fits.Column(name="DF",format='D', array=delta )
	hducols = md.columns + c0
	# now writes the catalog
	hdu = fits.BinTableHDU.from_columns(hducols)
	os.system("rm -rf "+path_to_outputCat)
	hdu.writeto(path_to_outputCat)

for el in snList:
	writeDFMock(el, DFfile, Lbox = 1000.)
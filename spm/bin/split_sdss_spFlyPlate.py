import os
import glob
import numpy as n

import astropy.io.fits as fits

files = os.path.join(os.environ['EBOSS_ROOT'], 'spectro', 'firefly', 'v1_0_1', '*', 'stellarpop', '*', 'spFlyPlate-*.fits')
#files = os.path.join(os.environ['EBOSS_ROOT'], 'spectro', 'firefly', 'v1_0_1', '26', 'stellarpop', '026*', 'spFlyPlate-*.fits')

splates = n.array(glob.glob(files))

def rewrite_plate_files(hd, mjd):
	path_2_out_file = os.path.join(os.path.dirname(splate), os.path.basename(splate)[:-5]+"-"+str(mjd)+".fits")
	selection = (hd[1].data['MJD'] == mjd)
	newtbdata = hd[1].data[selection]
	hdu_1 = fits.BinTableHDU(data=newtbdata)
		
	thdulist = fits.HDUList([hd[0], hdu_1])
	if os.path.isfile(path_2_out_file ):
		return 0. #os.remove(path_2_out_file )
	else :
		thdulist.writeto( path_2_out_file )
		return 1.
	
for splate in splates:
	hd = fits.open(splate)
	mjds = n.array(list(set(hd[1].data['MJD'])))
	for mjd in mjds :
		print splate, mjd
		rewrite_plate_files(hd, mjd)

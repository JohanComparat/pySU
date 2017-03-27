
from os.path import join
import os
import numpy as n
import glob 
import sys 
import time
import astropy.io.fits as fits


def concatenate_spFlyAll(root_dir, imf, dir, out_cat):
	spFlyCats = n.array( glob.glob( join( os.environ[root_dir], dir, "flyAll_catalogs","*.fits") ) )
	spFlyCats.sort()
	t0 = time.time()
	init_cat = spFlyCats[0]
	hdu_orig_table = fits.open(init_cat)
	orig_table = hdu_orig_table[1].data
	orig_cols = orig_table.columns
	table_all = orig_table

	for fitFile in spFlyCats:
		print fitFile
		nDATA = fits.open(fitFile)[1].data
		print nDATA.dtype
		print nDATA.shape
		try :#if len(fits.open(fitFile)[1].data)>0:
			table_all = n.hstack((table_all, nDATA))
		except TypeError:
			print "type error"
		
	newDat = n.transpose(table_all)
	new_cols = fits.ColDefs(newDat)
	hdu = fits.BinTableHDU.from_columns(new_cols)
	if os.path.isfile(out_cat):
		os.remove(out_cat)
	hdu.writeto(out_cat)
	print 0, time.time()-t0


imf='kr'
dir ='stellarpop-m11-kroupa'
out_cat = join( os.environ['SDSSDR12_DIR'], "catalogs", "FireflyGalaxyKroupaSdss26.fits")
#concatenate_spFlyAll('SDSSDR12_DIR', imf, dir , out_cat )

imf='ss'
dir ='stellarpop-m11-salpeter'
out_cat = join( os.environ['SDSSDR12_DIR'], "catalogs","FireflyGalaxysalpeterSdss26.fits")
concatenate_spFlyAll('SDSSDR12_DIR', imf, dir , out_cat )


imf='kr'
dir ='stellarpop-m11-kroupa-nodust'
out_cat = join( os.environ['SDSSDR12_DIR'], "catalogs","FireflyGalaxyKroupaNodustSdss26.fits")
concatenate_spFlyAll('SDSSDR12_DIR', imf, dir , out_cat )

imf='ss'
dir ='stellarpop-m11-salpeter-nodust'
out_cat = join( os.environ['SDSSDR12_DIR'], "catalogs","FireflyGalaxySalpeterNodustSdss26.fits")
concatenate_spFlyAll('SDSSDR12_DIR', imf, dir , out_cat )


imf='kr'
dir ='stellarpop-m11-kroupa'
out_cat = join( os.environ['EBOSSDR14_DIR'], "catalogs","FireflyGalaxyKroupaEbossDR14.fits")
concatenate_spFlyAll('EBOSSDR14_DIR', imf, dir , out_cat )

imf='ss'
dir ='stellarpop-m11-salpeter'
out_cat = join( os.environ['EBOSSDR14_DIR'], "catalogs","FireflyGalaxySalpeterEbossDR14.fits")
concatenate_spFlyAll('EBOSSDR14_DIR', imf, dir , out_cat )


imf='kr'
dir ='stellarpop-m11-kroupa-nodust'
out_cat = join( os.environ['EBOSSDR14_DIR'], "catalogs","FireflyGalaxyKroupaNoDustEbossDR14.fits")
concatenate_spFlyAll('EBOSSDR14_DIR', imf, dir , out_cat )

imf='ss'
dir ='stellarpop-m11-salpeter-nodust'
out_cat = join( os.environ['EBOSSDR14_DIR'], "catalogs","FireflyGalaxySalpeterNodustEbossDR14.fits")
concatenate_spFlyAll('EBOSSDR14_DIR', imf, dir , out_cat )


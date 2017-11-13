#! /usr/bin/env python
import sys
from os.path import join
import os
import time
import numpy as np
import glob

# for one galaxy spectrum
import GalaxySpectrumFIREFLY as gs

import astropy.io.fits as fits

catalog=fits.open(join(os.environ['DEEP2_DIR'], "catalogs", "zcat.deep2.dr4.v4.LFcatalogTC.Planck15.fits"))[1].data

out_file = join( os.environ['DEEP2_DIR'], 'catalogs', "zcat.deep2.dr4.v4.LFcatalogTC.Planck15_SNR.fits")

SNR = np.zeros(len(catalog))

def runSpec(catalog_entry, mask, objno):
	print catalog_entry['OBJNO'], catalog_entry['MASK'], catalog_entry['ZBEST'], catalog_entry['RA'],  catalog_entry['DEC']
	t0=time.time()
	path_to_spectrum = glob.glob(join(os.environ['DEEP2_DIR'], 'spectra', mask, '*', '*' + objno + '*_fc_tc.dat'))
	print path_to_spectrum
	if len(path_to_spectrum)>=1:
		spec=gs.GalaxySpectrumFIREFLY("-", milky_way_reddening=True)
		spec.openObservedDEEP2pectrum(catalog_entry)
		ok =(spec.bad_flags == 1) & (spec.flux>0) 
		return np.median(spec.flux[ok]/spec.error[ok])
	else:
		return -1.


for ii, catalog_entry in enumerate(catalog):
	mask=str(catalog_entry['MASK'])
	objno=str(catalog_entry['OBJNO'])
	SNR[ii]=runSpec(catalog_entry, mask, objno)
	print('SNR=',SNR[ii])
		

orig_cols = catalog.columns
all_cols = [fits.Column(name="SNR_MEDIAN_ALL", format='D', array=SNR)]
new_cols = fits.ColDefs(all_cols)

tbhdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)

prihdr = fits.Header()

prihdr['author'] = "JC"

prihdu = fits.PrimaryHDU(header=prihdr)

hdu = fits.HDUList([prihdu, tbhdu])

if os.path.isfile(out_file):
    os.remove(out_file)

hdu.writeto(out_file)


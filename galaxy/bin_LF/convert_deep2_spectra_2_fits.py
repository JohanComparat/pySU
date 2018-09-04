#! /usr/bin/env python

dat_2_month = {
  '01': 'jan',
  '02': 'feb',
  '03': 'mar',
  '04': 'apr',
  '05': 'may',
  '06': 'jui',
  '07': 'jul',
  '08': 'aug',
  '09': 'sep',
  '10': 'oct',
  '11': 'nov',
  '12': 'dec'
  }

import sys
from os.path import join
import os
import time
import numpy as np
import glob
import astropy.cosmology as co
cosmo = co.Planck13

# for one galaxy spectrum
import GalaxySpectrumFIREFLY as gs
import StellarPopulationModel as spm

import astropy.io.fits as fits

deep2_dir = '/home/comparat/data2/firefly/v1_1_0/DEEP2'
outputFolder = join( deep2_dir, 'raw_data', 'fits_spec')

catalog=fits.open(join( deep2_dir, "catalogs", "inputs/zcat.deep2.dr4.v4.fits.gz"))[1].data

catalog_entry = catalog[0]
date_split = catalog_entry['DATE'].split('-')
path_date = "".join(np.array([date_split[0], dat_2_month[date_split[1]] ,date_split[2]]))
path_folder = join( deep2_dir, 'raw_data/spectra', str(catalog_entry['MASK']), path_date)
path_to_spectrum = join(path_folder, 'spec1d.'+str(catalog_entry['MASK'])+"."+str(catalog_entry['SLIT']).zfill(3)+"."+str(catalog_entry['OBJNO']) )

def convert_spec_2_fits(catalog_entry, output_file, mask, objno):
        path_to_spectrum = glob.glob(join(deep2_dir, 'spectra', mask, '*', '*' + objno + '*_fc_tc.dat'))
        print path_to_spectrum
        
        if len(path_to_spectrum)>=1:
                spec=gs.GalaxySpectrumFIREFLY("-", milky_way_reddening=True)
                spec.openObservedDEEP2pectrum(catalog_entry)

                prihdr = fits.Header()
                prihdr['FILE']          = os.path.basename(output_file)
                prihdr['MASK']          = catalog_entry['MASK'] 
                prihdr['OBJNO']         = catalog_entry['OBJNO']   
                prihdr['RA']                    = catalog_entry['RA']
                prihdr['DEC']           = catalog_entry['DEC']
                prihdr['redshift']          = catalog_entry['ZBEST']
                prihdu = fits.PrimaryHDU(header=prihdr)

                waveCol = fits.Column(name="wavelength",format="D", unit="Angstrom", array= spec.wavelength)
                dataCol = fits.Column(name="flux",format="D", unit="1e-17erg/s/cm2/Angstrom", array= spec.flux)
                errorCol = fits.Column(name="flux_error",format="D", unit="1e-17erg/s/cm2/Angstrom", array= spec.error)
                
                cols = fits.ColDefs([ waveCol, dataCol, errorCol]) 
                tbhdu = fits.BinTableHDU.from_columns(cols)


                complete_hdus = fits.HDUList([prihdu, tbhdu])
                if os.path.isfile(output_file):
                        os.remove(output_file)
                complete_hdus.writeto(output_file)
        

print len(catalog), "N lines in the catalog"
for catalog_entry in catalog:
        mask=str(catalog_entry['MASK'])
        objno=str(catalog_entry['OBJNO'])
        output_file = join(outputFolder, 'deep2-'+mask+'-'+objno +".fits")
        #if os.path.isfile(output_file):
        #       print "pass", output_file
        #else:
        convert_spec_2_fits(catalog_entry, output_file, mask, objno)


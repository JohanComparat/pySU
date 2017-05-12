from os.path import join
import os
import numpy as n
import glob 
import sys 
import time
import astropy.io.fits as fits
from os.path import join
import astropy.cosmology as co
cosmo = co.Planck13
import astropy.io.fits as fits

# for one galaxy spectrum
import GalaxySpectrumFIREFLY as gs
import StellarPopulationModel as spm


env = 'DATA_DIR'
spec_dir = join( os.environ[env], "spm", "GAMAmock")
all_filenames = n.array(glob.glob(os.path.join(spec_dir, "gal_*.dat")))
all_filenames.sort()
print "N files=",len(all_filenames)

stellarpop_dir = join( os.environ[env], "spm", "GAMAmock", 'stellarpop')
out_dir = join( os.environ[env], "spm", "GAMAmock", 'results')
im_dir = os.path.join(os.environ[env], "spm", "GAMAmock", 'images')

imfs = n.array(["Chabrier", "Chabrier","Chabrier","Salpeter","Salpeter","Salpeter","Kroupa","Kroupa","Kroupa"])
libs = n.array(["miles", "elodie","stelib","miles","elodie","stelib","miles","elodie","stelib"])

path_2_out = lambda filename : os.path.join(out_dir, os.path.basename(filename)[:-4]+".fits")
path_2_im = lambda filename : os.path.join(im_dir, os.path.basename(filename)[:-4]+".png")

filenames = n.array(glob.glob(os.path.join(out_dir, "gal_*.fits")))
summary_catalog = join(spec_dir, "catalogs", "GAMA.mock.spectra.spm.fits")

headers_base =" ID" 
header_i = n.array(["age_lightW", "age_lightW_err_plus", "age_lightW_err_minus", "metallicity_lightW", "metallicity_lightW_err_plus", "metallicity_lightW_err_minus", "stellar_mass", "stellar_mass_err_plus", "stellar_mass_err_minus", "spm_EBV", "nComponentsSSP", "chi2red"])


dV=-9999.99
def get_table_entry_full(hduSPM):
	return n.array([ 10**hduSPM.header['age_lightW'], 10**hduSPM.header['age_lightW_up']-10**hduSPM.header['age_lightW'], 10**hduSPM.header['age_lightW']-10**hduSPM.header['age_lightW_low'], hduSPM.header['metallicity_lightW'], hduSPM.header['metallicity_lightW_up'] - hduSPM.header['metallicity_lightW'], hduSPM.header['metallicity_lightW'] - hduSPM.header['metallicity_lightW_low'], hduSPM.header['stellar_mass'], hduSPM.header['stellar_mass_up'] - hduSPM.header['stellar_mass'], hduSPM.header['stellar_mass'] - hduSPM.header['stellar_mass_low'], hduSPM.header['EBV'], hduSPM.header['ssp_number'], hduSPM.header['chi2']])
	

table_all = []
for filename in filenames:
	identifier = int(os.path.basename(filename).split('_')[1])
	hdus = fits.open(filename)
	table_entry = [[identifier]]
	header_final = [[" ID"]]
	for ii,hdu in enumerate(hdus[1:]):
		table_entry.append(get_table_entry_full( hdu ))
		add_me = n.array(["_"+imfs[ii]+"_"+libs[ii] for el in header_i])
		header_final.append(n.core.defchararray.add(header_i, add_me))
		#print header_final

	table_all.append(n.hstack((table_entry)))
	header = n.hstack((header_final))

newDat = n.transpose(table_all)

all_cols = []
for data_array, head in zip(newDat, header):
    all_cols.append(fits.Column(name=head, format='D', array=data_array))

new_cols = fits.ColDefs(all_cols)
hdu = fits.BinTableHDU.from_columns(new_cols)
if os.path.isfile(summary_catalog):
	os.remove(summary_catalog)

hdu.writeto(summary_catalog)


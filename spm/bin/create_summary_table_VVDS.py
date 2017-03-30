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

init_cat=join(os.environ['VVDS_DIR'], "catalogs", "VVDS_WIDE_summary.v1.fits")
plate_catalog = join(os.environ['VVDS_DIR'], "catalogs", "VVDS_WIDE_summary.v1.spm.fits")
hdu_orig_table = fits.open(init_cat)
orig_table = hdu_orig_table[1].data
orig_cols = orig_table.columns

kroupaFolder = join( os.environ['VVDS_DIR'], 'stellarpop-m11-kroupa', 'stellarpop')
salpeterFolder = join( os.environ['VVDS_DIR'], 'stellarpop-m11-salpeter', 'stellarpop')


dV=-9999.99
def get_table_entry_full(hduSPM):
	return n.array([ 10**hduSPM.header['age_lightW_mean'], 10**hduSPM.header['age_lightW_mean_up']-10**hduSPM.header['age_lightW_mean'], 10**hduSPM.header['age_lightW_mean']-10**hduSPM.header['age_lightW_mean_low'], hduSPM.header['metallicity_lightW_mean'], hduSPM.header['metallicity_lightW_mean_up'] - hduSPM.header['metallicity_lightW_mean'], hduSPM.header['metallicity_lightW_mean'] - hduSPM.header['metallicity_lightW_mean_low'], hduSPM.header['stellar_mass_mean'], hduSPM.header['stellar_mass_mean_up'] - hduSPM.header['stellar_mass_mean'], hduSPM.header['stellar_mass_mean'] - hduSPM.header['stellar_mass_mean_low'], hduSPM.header['EBV'], hduSPM.header['ssp_number']])
	
headers =" age_lightW_mean_kroupa age_lightW_err_plus_kroupa age_lightW_err_minus_kroupa metallicity_lightW_mean_kroupa metallicity_lightW_mean_err_plus_kroupa metallicity_lightW_mean_err_minus_kroupa stellar_mass_kroupa stellar_mass_err_plus_kroupa stellar_mass_err_minus_kroupa spm_EBV_kroupa nComponentsSSP_kroupa age_lightW_mean_salpeter age_lightW_err_plus_salpeter age_lightW_err_minus_salpeter metallicity_lightW_mean_salpeter metallicity_lightW_mean_err_plus_salpeter metallicity_lightW_mean_err_minus_salpeter stellar_mass_salpeter stellar_mass_err_plus_salpeter stellar_mass_err_minus_salpeter spm_EBV_salpeter nComponentsSSP_salpeter"

table_all = []
for catalog_entry in orig_table:
	krF = join(kroupaFolder, 'spFly-vvdswide-' + str(catalog_entry['NUM']) + "-kr.fits")
	ssF = join(salpeterFolder, 'spFly-vvdswide-' + str(catalog_entry['NUM']) + "-ss.fits")
	if os.path.isfile(krF) and os.path.isfile(ssF):
		#print "gets info"
		table_entry_kr = get_table_entry_full( hduSPM=fits.open(krF)[1] )
		#print table_entry_kr.shape
		table_entry_ss = get_table_entry_full( hduSPM=fits.open(ssF)[1] )
		#print table_entry_ss.shape
		table_entry = n.hstack((table_entry_kr, table_entry_ss))
		table_all.append(table_entry)
	else:
		table_all.append(n.ones(22)*dV)
		
newDat = n.transpose(table_all)

all_cols = []
for data_array, head in zip(newDat, headers.split()):
	all_cols.append(fits.Column(name=head, format='D', array=data_array))

new_cols = fits.ColDefs(all_cols)
hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
if os.path.isfile(plate_catalog):
	os.remove(plate_catalog)

hdu.writeto(plate_catalog)


init_cat=join(os.environ['VVDS_DIR'], "catalogs", "VVDS_DEEP_summary.v1.fits")
summary_catalog = join(os.environ['VVDS_DIR'], "catalogs", "VVDS_DEEP_summary.v1.spm.fits")
hdu_orig_table = fits.open(init_cat)
orig_table = hdu_orig_table[1].data
orig_cols = orig_table.columns


table_all = []
for catalog_entry in orig_table:
	krF = join(kroupaFolder,'spFly-vvdsdeep-'+str(catalog_entry['NUM'])+"-kr.fits")
	ssF = join(salpeterFolder, 'spFly-vvdsdeep-'+str(catalog_entry['NUM'])+"-ss.fits")
	if os.path.isfile(krF) and os.path.isfile(ssF):
		#print "gets info"
		table_entry_kr = get_table_entry_full( hduSPM=fits.open(krF)[1] )
		#print table_entry_kr.shape
		table_entry_ss = get_table_entry_full( hduSPM=fits.open(ssF)[1] )
		#print table_entry_ss.shape
		table_entry = n.hstack((table_entry_kr, table_entry_ss))
		table_all.append(table_entry)
	else:
		table_all.append(n.ones(22)*dV)
		
newDat = n.transpose(table_all)

all_cols = []
for data_array, head in zip(newDat, headers.split()):
	all_cols.append(fits.Column(name=head, format='D', array=data_array))

new_cols = fits.ColDefs(all_cols)
hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
if os.path.isfile(summary_catalog):
	os.remove(summary_catalog)

hdu.writeto(summary_catalog)

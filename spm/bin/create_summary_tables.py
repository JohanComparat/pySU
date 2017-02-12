
from os.path import join
import os
import numpy as n
import glob 
import sys 
import time
import astropy.io.fits as fits

plate = sys.argv[1]
dir ='stellarpop-m11-kroupa'

hdus = fits.open( join( os.environ['SDSSDR12_DIR'], "catalogs", "specObj-dr12.fits") )

#plates_all = n.loadtxt( join(os.environ['SDSSDR12_DIR'], "plateNumberList"), unpack=True, dtype='str')
#plates = plates_all[:-2]


dV=-9999.99
def get_table_entry_full(hduSPM):
	headerA =" age_universe age_lightW_mean age_lightW_err_plus age_lightW_err_minus metallicity_lightW_mean metallicity_lightW_mean_err_plus metallicity_lightW_mean_err_minus age_massW_mean age_massW_err_plus age_massW_err_minus metallicity_massW_mean metallicity_massW_mean_err_plus metallicity_massW_mean_err_minus stellar_mass stellar_mass_err_plus stellar_mass_err_minus spm_EBV nComponentsSSP"
	
	table_entry = [ 10**hduSPM.header['age_universe'], 10**hduSPM.header['age_lightW_mean'], 10**hduSPM.header['age_lightW_mean_up']-10**hduSPM.header['age_lightW_mean'], 10**hduSPM.header['age_lightW_mean']-10**hduSPM.header['age_lightW_mean_low'], hduSPM.header['metallicity_lightW_mean'], hduSPM.header['metallicity_lightW_mean_up'] - hduSPM.header['metallicity_lightW_mean'], hduSPM.header['metallicity_lightW_mean'] - hduSPM.header['metallicity_lightW_mean_low'], 10**hduSPM.header['age_massW_mean'], 10**hduSPM.header['age_massW_mean_up']-10**hduSPM.header['age_massW_mean'], 10**hduSPM.header['age_massW_mean']-10**hduSPM.header['age_massW_mean_low'], hduSPM.header['metallicity_massW_mean'], hduSPM.header['metallicity_massW_mean_up'] - hduSPM.header['metallicity_massW_mean'], hduSPM.header['metallicity_massW_mean'] - hduSPM.header['metallicity_massW_mean_low'], hduSPM.header['stellar_mass_mean'], hduSPM.header['stellar_mass_mean_up'] - hduSPM.header['stellar_mass_mean'], hduSPM.header['stellar_mass_mean'] - hduSPM.header['stellar_mass_mean_low'], hduSPM.header['EBV'], hduSPM.header['ssp_number']]
	#print hduSPM.header
	for iii in n.arange(hduSPM.header['ssp_number']):
		table_entry.append( hduSPM.header['stellar_mass_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['age_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['metal_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['SFR_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['weightMass_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['weightLight_ssp_'+str(iii)] )
		headerA += ' stellar_mass_ssp_'+str(iii) + ' age_ssp_'+str(iii) + ' metal_ssp_'+str(iii) + ' SFR_ssp_'+str(iii) + ' weightMass_ssp_'+str(iii) + ' weightLight_ssp_'+str(iii)
	
	if hduSPM.header['ssp_number']<8 :
		for iii in n.arange(hduSPM.header['ssp_number'], 8, 1):
			table_entry.append([dV, dV, dV, dV, dV, dV])
			headerA += ' stellar_mass_ssp_'+str(iii) + ' age_ssp_'+str(iii) + ' metal_ssp_'+str(iii) + ' SFR_ssp_'+str(iii) + ' weightMass_ssp_'+str(iii) + ' weightLight_ssp_'+str(iii)

	table_entry = n.array( n.hstack((table_entry)) )
	#print table_entry.shape
	return n.hstack((table_entry)), headerA
	
# step 2 : match to thecreated data set
init_cat = join( os.environ['SDSSDR12_DIR'], "catalogs", "perPlate", "sp-"+plate+".fits")
plate_catalog = join( os.environ['SDSSDR12_DIR'], dir, "catalogs", "spFly-"+plate+".fits")
	
hdu_orig_table = fits.open(init_cat)
orig_table = hdu_orig_table[1].data
orig_cols = orig_table.columns

table_all = []

for fiber, mjd in zip(orig_table['FIBERID'], orig_table['MJD']):
	fitFile = join( os.environ['SDSSDR12_DIR'], dir, "stellarpop", plate, "spec-"+plate+"-"+str(mjd)+"-"+str(fiber).zfill(4)+"-SPM-MILES.fits")
	if os.path.isfile(fitFile):
		table_entry, headers = get_table_entry_full( hduSPM=fits.open(fitFile)[1] )
		table_all.append(table_entry)
	else:
		table_all.append(n.ones(66)*dV)

newDat = n.transpose(table_all)

all_cols = []
for data_array, head in zip(newDat, headers.split()):
	all_cols.append(fits.Column(name=head, format='D', array=data_array))

new_cols = fits.ColDefs(all_cols)
hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
if os.path.isfile(plate_catalog):
	os.remove(plate_catalog)

hdu.writeto(plate_catalog)

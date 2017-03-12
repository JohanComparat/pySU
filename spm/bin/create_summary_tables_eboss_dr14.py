from os.path import join
import os
import numpy as n
import glob 
import sys 
import time
import astropy.io.fits as fits

plate = sys.argv[1] # '9003'
dir = sys.argv[2] # 'stellarpop-m11-kroupa'
#python create_summary_tables.py 2578 stellarpop-m11-kroupa
#python test.py 2578 stellarpop-m11-kroupa

if dir == 'stellarpop-m11-salpeter' or dir == 'stellarpop-m11-salpeter-nodust' :
	suffix = "-ss.fits"

if dir == 'stellarpop-m11-kroupa' or dir == 'stellarpop-m11-kroupa-nodust' :
	suffix = "-kr.fits"

print plate
# print dir
# print suffix

init_cat = join( os.environ['EBOSSDR14_DIR'], "catalogs", "perPlate", "sp-"+plate.zfill(4)+".fits")
plate_catalog = join( os.environ['EBOSSDR14_DIR'], dir, "catalogs", "spFly-"+plate.zfill(4)+".fits")

dV=-9999.99
def get_table_entry_full(hduSPM):
	# print "gets entry"
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
hdu_orig_table = fits.open(init_cat)
orig_table = hdu_orig_table[1].data
orig_cols = orig_table.columns

aaa=orig_cols.del_col('CHUNK'                      )
aaa=orig_cols.del_col('PLATEQUALITY'           )
aaa=orig_cols.del_col('PLATESN2'                  )
aaa=orig_cols.del_col('DEREDSN2'                )
aaa=orig_cols.del_col('LAMBDA_EFF'             )
aaa=orig_cols.del_col('BLUEFIBER'                )
aaa=orig_cols.del_col('ZOFFSET'                   )
aaa=orig_cols.del_col('SPECPRIMARY'            )
aaa=orig_cols.del_col('SPECBOSS'                 )
aaa=orig_cols.del_col('BOSS_SPECOBJ_ID'    )

aaa=orig_cols.del_col('NSPECOBS'                 )
aaa=orig_cols.del_col('CX'                             )
aaa=orig_cols.del_col('CY'                             )
aaa=orig_cols.del_col('CZ'                             )
aaa=orig_cols.del_col('XFOCAL'                     )
aaa=orig_cols.del_col('YFOCAL'                     )

aaa=orig_cols.del_col('TFILE'                         )
aaa=orig_cols.del_col('TCOLUMN'                  )
aaa=orig_cols.del_col('NPOLY'                       )
aaa=orig_cols.del_col('THETA'                       )

aaa=orig_cols.del_col('WAVEMIN'                  )
aaa=orig_cols.del_col('WAVEMAX'                  )
aaa=orig_cols.del_col('WCOVERAGE'             )

aaa=orig_cols.del_col('SN_MEDIAN_ALL'       )
aaa=orig_cols.del_col('SN_MEDIAN'               )
aaa=orig_cols.del_col('CHI68P'                      )
aaa=orig_cols.del_col('FRACNSIGMA'             )
aaa=orig_cols.del_col('FRACNSIGHI'              )
aaa=orig_cols.del_col('FRACNSIGLO'             )
aaa=orig_cols.del_col('SPECTROFLUX'           )
aaa=orig_cols.del_col('SPECTROFLUX_IVAR'  )
aaa=orig_cols.del_col('SPECTROSYNFLUX'     )
aaa=orig_cols.del_col('SPECTROSYNFLUX_IVAR'  )
aaa=orig_cols.del_col('SPECTROSKYFLUX'           )
aaa=orig_cols.del_col('ANYANDMASK'                  )
aaa=orig_cols.del_col('ANYORMASK'                    )
aaa=orig_cols.del_col('SPEC1_G'                         )
aaa=orig_cols.del_col('SPEC1_R'                          )
aaa=orig_cols.del_col('SPEC1_I'                          )
aaa=orig_cols.del_col('SPEC2_G'                         )
aaa=orig_cols.del_col('SPEC2_R'                         ) 
aaa=orig_cols.del_col('SPEC2_I'                          )
aaa=orig_cols.del_col('ELODIE_FILENAME'          )
aaa=orig_cols.del_col('ELODIE_OBJECT'              )
aaa=orig_cols.del_col('ELODIE_SPTYPE'               )
aaa=orig_cols.del_col('ELODIE_BV'                      )
aaa=orig_cols.del_col('ELODIE_TEFF'                   )
aaa=orig_cols.del_col('ELODIE_LOGG'                 )
aaa=orig_cols.del_col('ELODIE_FEH'                    )
aaa=orig_cols.del_col('ELODIE_Z'                        )
aaa=orig_cols.del_col('ELODIE_Z_ERR'                )
aaa=orig_cols.del_col('ELODIE_Z_MODELERR'     )
aaa=orig_cols.del_col('ELODIE_RCHI2'                )
aaa=orig_cols.del_col('ELODIE_DOF'                   )

aaa=orig_cols.del_col('CALIBFLUX'                      )
aaa=orig_cols.del_col('CALIBFLUX_IVAR'             )

table_all = []
headers = ""
for fiber, mjd in zip(orig_table['FIBERID'], orig_table['MJD']):
	fitFile = join( os.environ['EBOSSDR14_DIR'], dir, "stellarpop", plate, "spFly-"+plate.zfill(4)+"-"+str(mjd)+"-"+str(fiber).zfill(4)+suffix)
	# print fitFile
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


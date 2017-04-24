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

init_cat = join( os.environ['SDSSDR12_DIR'], "catalogs", "perPlate", "sp-"+plate.zfill(4)+".fits")
plate_catalog = join( os.environ['SDSSDR12_DIR'], dir, "catalogs", "spFly-"+plate.zfill(4)+".fits")

dV=-9999.99
def get_table_entry_full(hduSPM):
	# print "gets entry"
	headerA =" age_universe age_lightW age_lightW_err_plus age_lightW_err_minus metallicity_lightW metallicity_lightW_err_plus metallicity_lightW_err_minus age_massW age_massW_err_plus age_massW_err_minus metallicity_massW metallicity_massW_err_plus metallicity_massW_err_minus stellar_mass stellar_mass_err_plus stellar_mass_err_minus spm_EBV nComponentsSSP"
	
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

orig_cols.del_col('CHUNK'                      )
#orig_cols.del_col('PROGRAMNAME'         )

orig_cols.del_col('PLATERUN'                 )
#orig_cols.del_col('PLATEQUALITY'           )
#orig_cols.del_col('PLATESN2'                  )
orig_cols.del_col('DEREDSN2'                )
orig_cols.del_col('LAMBDA_EFF'             )
orig_cols.del_col('BLUEFIBER'                )
orig_cols.del_col('ZOFFSET'                   )
orig_cols.del_col('SNTURNOFF'              )
orig_cols.del_col('NTURNOFF'                )
orig_cols.del_col('SPECPRIMARY'            )
orig_cols.del_col('SPECSDSS'                 )
orig_cols.del_col('SPECLEGACY'             )
orig_cols.del_col('SPECSEGUE'               )
orig_cols.del_col('SPECSEGUE1'             )
orig_cols.del_col('SPECSEGUE2'             )
orig_cols.del_col('SPECBOSS'                 )
#orig_cols.del_col('BOSS_SPECOBJ_ID'    )

orig_cols.del_col('SEGUE1_TARGET1'     )
orig_cols.del_col('SEGUE1_TARGET2'     )
orig_cols.del_col('SEGUE2_TARGET1'     )
orig_cols.del_col('SEGUE2_TARGET2'     )
orig_cols.del_col('MARVELS_TARGET1'   )
orig_cols.del_col('MARVELS_TARGET2'   )

orig_cols.del_col('PLATEID'                    )
orig_cols.del_col('NSPECOBS'                 )
orig_cols.del_col('FIRSTRELEASE'           )
orig_cols.del_col('DESIGNID'                 )
orig_cols.del_col('CX'                             )
orig_cols.del_col('CY'                             )
orig_cols.del_col('CZ'                             )
orig_cols.del_col('XFOCAL'                     )
orig_cols.del_col('YFOCAL'                     )

orig_cols.del_col('TFILE'                         )
orig_cols.del_col('TCOLUMN'                  )
orig_cols.del_col('NPOLY'                       )
orig_cols.del_col('THETA'                       )

orig_cols.del_col('WAVEMIN'                  )
orig_cols.del_col('WAVEMAX'                  )
orig_cols.del_col('WCOVERAGE'             )

#orig_cols.del_col('SN_MEDIAN_ALL'       )
#orig_cols.del_col('SN_MEDIAN'               )
orig_cols.del_col('CHI68P'                      )
orig_cols.del_col('FRACNSIGMA'             )
orig_cols.del_col('FRACNSIGHI'              )
orig_cols.del_col('FRACNSIGLO'             )
orig_cols.del_col('SPECTROFLUX'           )
orig_cols.del_col('SPECTROFLUX_IVAR'  )
orig_cols.del_col('SPECTROSYNFLUX'     )
orig_cols.del_col('SPECTROSYNFLUX_IVAR'  )
orig_cols.del_col('SPECTROSKYFLUX'           )
orig_cols.del_col('ANYANDMASK'                  )
orig_cols.del_col('ANYORMASK'                    )
orig_cols.del_col('SPEC1_G'                         )
orig_cols.del_col('SPEC1_R'                          )
orig_cols.del_col('SPEC1_I'                          )
orig_cols.del_col('SPEC2_G'                         )
orig_cols.del_col('SPEC2_R'                         ) 
orig_cols.del_col('SPEC2_I'                          )
orig_cols.del_col('ELODIE_FILENAME'          )
orig_cols.del_col('ELODIE_OBJECT'              )
orig_cols.del_col('ELODIE_SPTYPE'               )
orig_cols.del_col('ELODIE_BV'                      )
orig_cols.del_col('ELODIE_TEFF'                   )
orig_cols.del_col('ELODIE_LOGG'                 )
orig_cols.del_col('ELODIE_FEH'                    )
orig_cols.del_col('ELODIE_Z'                        )
orig_cols.del_col('ELODIE_Z_ERR'                )
orig_cols.del_col('ELODIE_Z_MODELERR'     )
orig_cols.del_col('ELODIE_RCHI2'                )
orig_cols.del_col('ELODIE_DOF'                   )

orig_cols.del_col('Z_PERSON'                       )
orig_cols.del_col('CLASS_PERSON'               )
orig_cols.del_col('Z_CONF_PERSON'            ) 
orig_cols.del_col('COMMENTS_PERSON'       )
orig_cols.del_col('CALIBFLUX'                      )
orig_cols.del_col('CALIBFLUX_IVAR'             )

table_all = []
headers = ""
for fiber, mjd in zip(orig_table['FIBERID'], orig_table['MJD']):
	fitFile = join( os.environ['SDSSDR12_DIR'], dir, "stellarpop", plate, "spFly-"+plate.zfill(4)+"-"+str(mjd)+"-"+str(fiber).zfill(4)+suffix)
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


sys.exit()



orig_cols.del_col('SURVEY'                     )
orig_cols.del_col('INSTRUMENT'             )
orig_cols.del_col('CHUNK'                      )
orig_cols.del_col('PROGRAMNAME'         )
orig_cols.del_col('PLATERUN'                 )
orig_cols.del_col('PLATEQUALITY'           )
orig_cols.del_col('PLATESN2'                  )
orig_cols.del_col('DEREDSN2'                )
orig_cols.del_col('LAMBDA_EFF'             )
orig_cols.del_col('BLUEFIBER'                )
orig_cols.del_col('ZOFFSET'                   )
orig_cols.del_col('SNTURNOFF'              )
orig_cols.del_col('NTURNOFF'                )
orig_cols.del_col('SPECPRIMARY'            )
orig_cols.del_col('SPECSDSS'                 )
orig_cols.del_col('SPECLEGACY'             )
orig_cols.del_col('SPECSEGUE'               )
orig_cols.del_col('SPECSEGUE1'             )
orig_cols.del_col('SPECSEGUE2'             )
orig_cols.del_col('SPECBOSS'                 )
orig_cols.del_col('BOSS_SPECOBJ_ID'    )
orig_cols.del_col('SPECOBJID'                )
orig_cols.del_col('FLUXOBJID'                )
orig_cols.del_col('BESTOBJID'                )
orig_cols.del_col('TARGETOBJID'           )
orig_cols.del_col('PLATEID'                    )
orig_cols.del_col('NSPECOBS'                 )
orig_cols.del_col('FIRSTRELEASE'           )
orig_cols.del_col('RUN2D'                      )
orig_cols.del_col('RUN1D'                      )
orig_cols.del_col('DESIGNID'                 )
orig_cols.del_col('CX'                             )
orig_cols.del_col('CY'                             )
orig_cols.del_col('CZ'                             )
orig_cols.del_col('XFOCAL'                     )
orig_cols.del_col('YFOCAL'                     )
orig_cols.del_col('SOURCETYPE'             )
orig_cols.del_col('TARGETTYPE'             )
orig_cols.del_col('PRIMTARGET'             )
orig_cols.del_col('SECTARGET'               )
orig_cols.del_col('LEGACY_TARGET1'     )
orig_cols.del_col('LEGACY_TARGET2'     )
orig_cols.del_col('SPECIAL_TARGET1'    )
orig_cols.del_col('SPECIAL_TARGET2'    )
orig_cols.del_col('SEGUE1_TARGET1'     )
orig_cols.del_col('SEGUE1_TARGET2'     )
orig_cols.del_col('SEGUE2_TARGET1'     )
orig_cols.del_col('SEGUE2_TARGET2'     )
orig_cols.del_col('MARVELS_TARGET1'   )
orig_cols.del_col('MARVELS_TARGET2'   )
orig_cols.del_col('BOSS_TARGET1'         )
orig_cols.del_col('BOSS_TARGET2'         )
orig_cols.del_col('EBOSS_TARGET0'       )
orig_cols.del_col('ANCILLARY_TARGET1')
orig_cols.del_col('ANCILLARY_TARGET2')
orig_cols.del_col('SPECTROGRAPHID'     )
orig_cols.del_col('PLATE'                        )
orig_cols.del_col('TILE'                          )
orig_cols.del_col('MJD'                           )
orig_cols.del_col('FIBERID'                    )
orig_cols.del_col('OBJID'                        )
orig_cols.del_col('PLUG_RA'                   )
orig_cols.del_col('PLUG_DEC'                 )
orig_cols.del_col('CLASS'                       )
orig_cols.del_col('SUBCLASS'                 )
orig_cols.del_col('Z'                               )	
orig_cols.del_col('Z_ERR'                       )
orig_cols.del_col('RCHI2'                       )
orig_cols.del_col('DOF'                          )
orig_cols.del_col('RCHI2DIFF'                )
orig_cols.del_col('TFILE'                         )
orig_cols.del_col('TCOLUMN'                  )
orig_cols.del_col('NPOLY'                       )
orig_cols.del_col('THETA'                       )
orig_cols.del_col('VDISP'                        )
orig_cols.del_col('VDISP_ERR'                )
orig_cols.del_col('VDISPZ'                      )
orig_cols.del_col('VDISPZ_ERR'              )
orig_cols.del_col('VDISPCHI2'                )
orig_cols.del_col('VDISPNPIX'                )
orig_cols.del_col('VDISPDOF'                 )
orig_cols.del_col('WAVEMIN'                  )
orig_cols.del_col('WAVEMAX'                  )
orig_cols.del_col('WCOVERAGE'             )
orig_cols.del_col('ZWARNING'                )
orig_cols.del_col('SN_MEDIAN_ALL'       )
orig_cols.del_col('SN_MEDIAN'               )
orig_cols.del_col('CHI68P'                      )
orig_cols.del_col('FRACNSIGMA'             )
orig_cols.del_col('FRACNSIGHI'              )
orig_cols.del_col('FRACNSIGLO'             )
orig_cols.del_col('SPECTROFLUX'           )
orig_cols.del_col('SPECTROFLUX_IVAR'  )
orig_cols.del_col('SPECTROSYNFLUX'     )
orig_cols.del_col('SPECTROSYNFLUX_IVAR'  )
orig_cols.del_col('SPECTROSKYFLUX'           )
orig_cols.del_col('ANYANDMASK'                  )
orig_cols.del_col('ANYORMASK'                    )
orig_cols.del_col('SPEC1_G'                         )
orig_cols.del_col('SPEC1_R'                          )
orig_cols.del_col('SPEC1_I'                          )
orig_cols.del_col('SPEC2_G'                         )
orig_cols.del_col('SPEC2_R'                         ) 
orig_cols.del_col('SPEC2_I'                          )
orig_cols.del_col('ELODIE_FILENAME'          )
orig_cols.del_col('ELODIE_OBJECT'              )
orig_cols.del_col('ELODIE_SPTYPE'               )
orig_cols.del_col('ELODIE_BV'                      )
orig_cols.del_col('ELODIE_TEFF'                   )
orig_cols.del_col('ELODIE_LOGG'                 )
orig_cols.del_col('ELODIE_FEH'                    )
orig_cols.del_col('ELODIE_Z'                        )
orig_cols.del_col('ELODIE_Z_ERR'                )
orig_cols.del_col('ELODIE_Z_MODELERR'     )
orig_cols.del_col('ELODIE_RCHI2'                )
orig_cols.del_col('ELODIE_DOF'                   )
orig_cols.del_col('Z_NOQSO'                       )
orig_cols.del_col('Z_ERR_NOQSO'               )
orig_cols.del_col('ZWARNING_NOQSO'        )
orig_cols.del_col('CLASS_NOQSO'                )
orig_cols.del_col('SUBCLASS_NOQSO'          )
orig_cols.del_col('RCHI2DIFF_NOQSO'         )
orig_cols.del_col('Z_PERSON'                       )
orig_cols.del_col('CLASS_PERSON'               )
orig_cols.del_col('Z_CONF_PERSON'            ) 
orig_cols.del_col('COMMENTS_PERSON'       )
orig_cols.del_col('CALIBFLUX'                      )
orig_cols.del_col('CALIBFLUX_IVAR'             )

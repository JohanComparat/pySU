
from os.path import join
import os
import numpy as n
import glob 
import sys 
import time
import astropy.io.fits as fits


def concatenate_spFlyPlates(root_dir, imf, dir, all_cat):
	spFlyCats = n.array( glob.glob( join( os.environ[root_dir], dir, "catalogs","spFly-*.fits") ) )
	spFlyCats.sort()

	NperCat = 100
	bds = n.arange(0, len(spFlyCats), NperCat)
	for bd in bds:
		t0 = time.time()
		init_cat = spFlyCats[bd]
		hdu_orig_table = fits.open(init_cat)
		orig_table = hdu_orig_table[1].data
		orig_cols = orig_table.columns

		table_all = orig_table

		for fitFile in spFlyCats[bd+1:bd+NperCat]:
			print fitFile
			nDATA = fits.open(fitFile)[1].data
			#print nDATA
			try :#if len(fits.open(fitFile)[1].data)>0:
				table_all = n.hstack((table_all, nDATA))
			except TypeError:
				print "type error"
			
		newDat = n.transpose(table_all)

		new_cols = fits.ColDefs(newDat)

		hdu = fits.BinTableHDU.from_columns(new_cols)
		write_cat = all_cat+"-"+str(bd)+".fits"
		if os.path.isfile(write_cat):
			os.remove(write_cat)

		hdu.writeto(write_cat)
		print bd, time.time()-t0

imf='kr'
dir ='stellarpop-m11-kroupa'
all_cat = join( os.environ['EBOSSDR14_DIR'], dir, "flyAll_catalogs", "spFlyAll-"+imf)
concatenate_spFlyPlates('EBOSSDR14_DIR', imf, dir , all_cat )

imf='ss'
dir ='stellarpop-m11-salpeter'
all_cat = join( os.environ['EBOSSDR14_DIR'], dir, "flyAll_catalogs", "spFlyAll-"+imf)
concatenate_spFlyPlates('EBOSSDR14_DIR', imf, dir , all_cat )


imf='kr'
dir ='stellarpop-m11-kroupa-nodust'
all_cat = join( os.environ['EBOSSDR14_DIR'], dir, "flyAll_catalogs", "spFlyAll-"+imf)
concatenate_spFlyPlates('EBOSSDR14_DIR', imf, dir , all_cat )

imf='ss'
dir ='stellarpop-m11-salpeter-nodust'
all_cat = join( os.environ['EBOSSDR14_DIR'], dir, "flyAll_catalogs", "spFlyAll-"+imf)
concatenate_spFlyPlates('EBOSSDR14_DIR', imf, dir , all_cat )


imf='kr'
dir ='stellarpop-m11-kroupa'
all_cat = join( os.environ['SDSSDR12_DIR'], dir, "flyAll_catalogs", "spFlyAll-"+imf)
concatenate_spFlyPlates('SDSSDR12_DIR', imf, dir , all_cat )

imf='ss'
dir ='stellarpop-m11-salpeter'
all_cat = join( os.environ['SDSSDR12_DIR'], dir, "flyAll_catalogs", "spFlyAll-"+imf)
concatenate_spFlyPlates('SDSSDR12_DIR', imf, dir , all_cat )


imf='kr'
dir ='stellarpop-m11-kroupa-nodust'
all_cat = join( os.environ['SDSSDR12_DIR'], dir, "flyAll_catalogs", "spFlyAll-"+imf)
concatenate_spFlyPlates('SDSSDR12_DIR', imf, dir , all_cat )

imf='ss'
dir ='stellarpop-m11-salpeter-nodust'
all_cat = join( os.environ['SDSSDR12_DIR'], dir, "flyAll_catalogs", "spFlyAll-"+imf)
concatenate_spFlyPlates('SDSSDR12_DIR', imf, dir , all_cat )


sys.exit()

headers = n.array(['specObjID','mjd','plate','fiberID','run1d','run2d','ra','dec','z_noqso','zErr_noqso','zWarning_noqso','class_noqso','subClass_noqso','u','g','r','i','z','err_u','err_g','err_r','err_i','err_z','dered_u','dered_g','dered_r','dered_i','dered_z','age_mean','age_err_plus','age_err_minus','metallicity_mean','metallicity_mean_err_plus','metallicity_mean_err_minus','stellar_mass','stellar_mass_err_plus	','stellar_mass_err_minus','spm_EBV','nComponentsSSP'])

new_cols.del_col('SURVEY'             )
new_cols.del_col('INSTRUMENT'     )
#new_cols.del_col('CHUNK'              )
new_cols.del_col('PROGRAMNAME' )
new_cols.del_col('PLATERUN'         )
new_cols.del_col('PLATEQUALITY'   )
new_cols.del_col('PLATESN2'          )
new_cols.del_col('DEREDSN2'        )
new_cols.del_col('LAMBDA_EFF'     )
new_cols.del_col('BLUEFIBER'        )
new_cols.del_col('ZOFFSET')
new_cols.del_col('SNTURNOFF'      )
new_cols.del_col('NTURNOFF'        )
new_cols.del_col('SPECPRIMARY'    )
new_cols.del_col('SPECSDSS'         )
new_cols.del_col('SPECLEGACY'     )
new_cols.del_col('SPECSEGUE'       )
new_cols.del_col('SPECSEGUE1'     )
new_cols.del_col('SPECSEGUE2'     )
new_cols.del_col('SPECBOSS'         )
new_cols.del_col('BOSS_SPECOBJ_ID'    )
#new_cols.del_col('SPECOBJID'                )
new_cols.del_col('FLUXOBJID'      )
new_cols.del_col('BESTOBJID'      )
new_cols.del_col('TARGETOBJID' )
new_cols.del_col('PLATEID'          )
new_cols.del_col('NSPECOBS'       )
new_cols.del_col('FIRSTRELEASE' )
#new_cols.del_col('RUN2D'            )
#new_cols.del_col('RUN1D'            )
new_cols.del_col('DESIGNID'       )
new_cols.del_col('CX'                   )
new_cols.del_col('CY'                   )
new_cols.del_col('CZ'                   )
new_cols.del_col('XFOCAL')
new_cols.del_col('YFOCAL')
new_cols.del_col('SOURCETYPE'   )
new_cols.del_col('TARGETTYPE'   )
new_cols.del_col('PRIMTARGET'   )
new_cols.del_col('SECTARGET'     )
new_cols.del_col('LEGACY_TARGET1'     )
new_cols.del_col('LEGACY_TARGET2'     )
new_cols.del_col('SPECIAL_TARGET1'    )
new_cols.del_col('SPECIAL_TARGET2'    )
new_cols.del_col('SEGUE1_TARGET1'     )
new_cols.del_col('SEGUE1_TARGET2'     )
new_cols.del_col('SEGUE2_TARGET1'     )
new_cols.del_col('SEGUE2_TARGET2'     )
new_cols.del_col('MARVELS_TARGET1'   )
new_cols.del_col('MARVELS_TARGET2'   )
new_cols.del_col('BOSS_TARGET1'         )
new_cols.del_col('BOSS_TARGET2'         )
new_cols.del_col('EBOSS_TARGET0'       )
new_cols.del_col('ANCILLARY_TARGET1')
new_cols.del_col('ANCILLARY_TARGET2')
new_cols.del_col('SPECTROGRAPHID'     )
#new_cols.del_col('PLATE'            )
new_cols.del_col('TILE')
#new_cols.del_col('MJD'               )
#new_cols.del_col('FIBERID'       )
new_cols.del_col('OBJID'            )
#new_cols.del_col('PLUG_RA'       )
#new_cols.del_col('PLUG_DEC'     )
new_cols.del_col('CLASS')
new_cols.del_col('SUBCLASS'    )
new_cols.del_col('Z')
new_cols.del_col('Z_ERR'          )
new_cols.del_col('RCHI2')
new_cols.del_col('DOF'              )
new_cols.del_col('RCHI2DIFF'    )
new_cols.del_col('TFILE'             )
new_cols.del_col('TCOLUMN'     )
new_cols.del_col('NPOLY')
new_cols.del_col('THETA')
new_cols.del_col('VDISP'            )
new_cols.del_col('VDISP_ERR'    )
new_cols.del_col('VDISPZ'          )
new_cols.del_col('VDISPZ_ERR'  )
new_cols.del_col('VDISPCHI2'    )
new_cols.del_col('VDISPNPIX'    )
new_cols.del_col('VDISPDOF'     )
new_cols.del_col('WAVEMIN'      )
new_cols.del_col('WAVEMAX'      )
new_cols.del_col('WCOVERAGE'             )
new_cols.del_col('ZWARNING'                )
new_cols.del_col('SN_MEDIAN_ALL'       )
new_cols.del_col('SN_MEDIAN'               )
new_cols.del_col('CHI68P'                      )
new_cols.del_col('FRACNSIGMA'             )
new_cols.del_col('FRACNSIGHI'              )
new_cols.del_col('FRACNSIGLO'             )
new_cols.del_col('SPECTROFLUX')
new_cols.del_col('SPECTROFLUX_IVAR'  )
new_cols.del_col('SPECTROSYNFLUX'     )
new_cols.del_col('SPECTROSYNFLUX_IVAR'  )
new_cols.del_col('SPECTROSKYFLUX')
new_cols.del_col('ANYANDMASK'                  )
new_cols.del_col('ANYORMASK'                    )
new_cols.del_col('SPEC1_G'                         )
new_cols.del_col('SPEC1_R'                          )
new_cols.del_col('SPEC1_I'                          )
new_cols.del_col('SPEC2_G'                         )
new_cols.del_col('SPEC2_R'                         ) 
new_cols.del_col('SPEC2_I'                          )
new_cols.del_col('ELODIE_FILENAME'          )
new_cols.del_col('ELODIE_OBJECT'              )
new_cols.del_col('ELODIE_SPTYPE'               )
new_cols.del_col('ELODIE_BV'                      )
new_cols.del_col('ELODIE_TEFF'                   )
new_cols.del_col('ELODIE_LOGG'                 )
new_cols.del_col('ELODIE_FEH'                    )
new_cols.del_col('ELODIE_Z'                        )
new_cols.del_col('ELODIE_Z_ERR'                )
new_cols.del_col('ELODIE_Z_MODELERR'     )
new_cols.del_col('ELODIE_RCHI2'                )
new_cols.del_col('ELODIE_DOF'                   )
#new_cols.del_col('Z_NOQSO'                       )
#new_cols.del_col('Z_ERR_NOQSO'               )
#new_cols.del_col('ZWARNING_NOQSO'        )
#new_cols.del_col('CLASS_NOQSO'                )
#new_cols.del_col('SUBCLASS_NOQSO'          )
new_cols.del_col('RCHI2DIFF_NOQSO'         )
new_cols.del_col('Z_PERSON'                       )
new_cols.del_col('CLASS_PERSON'               )
new_cols.del_col('Z_CONF_PERSON'            ) 
new_cols.del_col('COMMENTS_PERSON'       )
new_cols.del_col('CALIBFLUX'                      )
new_cols.del_col('CALIBFLUX_IVAR'             )


























orig_cols.del_col('CHUNK'                      )
orig_cols.del_col('PROGRAMNAME'         )

orig_cols.del_col('PLATERUN'                 )
orig_cols.del_col('PLATEQUALITY')
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

orig_cols.del_col('SEGUE1_TARGET1'     )
orig_cols.del_col('SEGUE1_TARGET2'     )
orig_cols.del_col('SEGUE2_TARGET1'     )
orig_cols.del_col('SEGUE2_TARGET2'     )
orig_cols.del_col('MARVELS_TARGET1'   )
orig_cols.del_col('MARVELS_TARGET2'   )

orig_cols.del_col('PLATEID'                    )
orig_cols.del_col('NSPECOBS'                 )
orig_cols.del_col('FIRSTRELEASE')
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

orig_cols.del_col('SN_MEDIAN_ALL'       )
orig_cols.del_col('SN_MEDIAN'               )
orig_cols.del_col('CHI68P'                      )
orig_cols.del_col('FRACNSIGMA'             )
orig_cols.del_col('FRACNSIGHI'              )
orig_cols.del_col('FRACNSIGLO'             )
orig_cols.del_col('SPECTROFLUX')
orig_cols.del_col('SPECTROFLUX_IVAR'  )
orig_cols.del_col('SPECTROSYNFLUX'     )
orig_cols.del_col('SPECTROSYNFLUX_IVAR'  )
orig_cols.del_col('SPECTROSKYFLUX')
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


orig_cols.del_col('SURVEY'                     )
orig_cols.del_col('INSTRUMENT'             )
orig_cols.del_col('CHUNK'                      )
orig_cols.del_col('PROGRAMNAME'         )
orig_cols.del_col('PLATERUN'                 )
orig_cols.del_col('PLATEQUALITY')
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
orig_cols.del_col('TARGETOBJID')
orig_cols.del_col('PLATEID'                    )
orig_cols.del_col('NSPECOBS'                 )
orig_cols.del_col('FIRSTRELEASE')
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
orig_cols.del_col('SPECTROFLUX')
orig_cols.del_col('SPECTROFLUX_IVAR'  )
orig_cols.del_col('SPECTROSYNFLUX'     )
orig_cols.del_col('SPECTROSYNFLUX_IVAR'  )
orig_cols.del_col('SPECTROSKYFLUX')
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

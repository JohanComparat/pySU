import time
t0t=time.time()
from os.path import join
import os
import numpy as n
import glob 
import sys 
import astropy.io.fits as fits

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

dV=-9999.99

plate = sys.argv[1] # '9003'
env = sys.argv[2] # '9003'
print plate, env
#python create_summary_tables_sdss_dr12.py 3785 EBOSSDR14_DIR

# initial catalog
init_cat = join( os.environ[env], "catalogs", "perPlate", "sp-"+plate.zfill(4)+".fits")

dir = 'stellarpop' 
suffix = ".fits"

out_dir = os.path.join(os.environ[env], 'stellarpop', plate)
im_dir = os.path.join(os.environ[env], 'stellarpop', plate, 'images')

if os.path.isdir(out_dir)==False:
	os.makedirs(out_dir)
if os.path.isdir(im_dir)==False:
	os.makedirs(im_dir)
	
path_2_out_file = join( out_dir, "spFlyPlate-"+plate.zfill(4)+".fits")

im_file = "spFlyPlate-"+plate.zfill(4)+".png"
path_2_im_file = os.path.join(im_dir, im_file)

path_2_im_file

prihdr = fits.Header()
prihdr['file']   = os.path.basename(path_2_out_file)
prihdr['plate']  = int(plate)
prihdr['models'] = 'Maraston_2011'
prihdr['library'] = 'MILES'
prihdr['fitter'] = 'FIREFLY'
prihdr['author'] = 'johan comparat'
prihdr['DR'] = 14


def get_table_entry_full(hduSPM):
	# print "gets entry"
	hduSPM.header
	prefix = hduSPM.header['IMF'] + "_" #+ hduSPM.header['library'] + "_"
	#print prefix
	headerA =" "+prefix+"age_lightW "+prefix+"age_lightW_err_plus "+prefix+"age_lightW_err_minus "+prefix+"metallicity_lightW "+prefix+"metallicity_lightW_err_plus "+prefix+"metallicity_lightW_err_minus "+prefix+"age_massW "+prefix+"age_massW_err_plus "+prefix+"age_massW_err_minus "+prefix+"metallicity_massW "+prefix+"metallicity_massW_err_plus "+prefix+"metallicity_massW_err_minus "+prefix+"stellar_mass "+prefix+"stellar_mass_err_plus "+prefix+"stellar_mass_err_minus "+prefix+"spm_EBV "+prefix+"nComponentsSSP"
	
	table_entry = [10**hduSPM.header['age_lightW']          
	, 10**hduSPM.header['age_lightW_up']       
	, 10**hduSPM.header['age_lightW_low']      
	, hduSPM.header['metallicity_lightW']  
	, hduSPM.header['metallicity_lightW_up']
	, hduSPM.header['metallicity_lightW_low']
	, 10**hduSPM.header['age_massW']           
	, 10**hduSPM.header['age_massW_up']        
	, 10**hduSPM.header['age_massW_low']       
	, hduSPM.header['metallicity_massW']   
	, hduSPM.header['metallicity_massW_up']
	, hduSPM.header['metallicity_massW_low']
	, hduSPM.header['stellar_mass']        
	, hduSPM.header['stellar_mass_up']     
	, hduSPM.header['stellar_mass_low']    
	, hduSPM.header['EBV'] 
	, hduSPM.header['ssp_number']
	]
	
	#print hduSPM.header
	for iii in n.arange(hduSPM.header['ssp_number']):
		table_entry.append( hduSPM.header['stellar_mass_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['age_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['metal_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['SFR_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['weightMass_ssp_'+str(iii)] )
		table_entry.append( hduSPM.header['weightLight_ssp_'+str(iii)] )
		headerA += ' '+prefix+'stellar_mass_ssp_'+str(iii) + ' '+prefix+'age_ssp_'+str(iii) + ' '+prefix+'metal_ssp_'+str(iii) + ' '+prefix+'SFR_ssp_'+str(iii) + ' '+prefix+'weightMass_ssp_'+str(iii) + ' '+prefix+'weightLight_ssp_'+str(iii)
	
	if hduSPM.header['ssp_number']<8 :
		for iii in n.arange(hduSPM.header['ssp_number'], 8, 1):
			table_entry.append([dV, dV, dV, dV, dV, dV])
			headerA += ' '+prefix+'stellar_mass_ssp_'+str(iii) + ' '+prefix+'age_ssp_'+str(iii) + ' '+prefix+'metal_ssp_'+str(iii) + ' '+prefix+'SFR_ssp_'+str(iii) + ' '+prefix+'weightMass_ssp_'+str(iii) + ' '+prefix+'weightLight_ssp_'+str(iii)

	table_entry = n.array( n.hstack((table_entry)) )
	#print table_entry.shape
	return n.hstack((table_entry)), headerA
	
# step 2 : match to thecreated data set	
hdu_orig_table = fits.open(init_cat)
orig_table = hdu_orig_table[1].data
orig_cols = orig_table.columns

table_all = []
headers = ""
for fiber, mjd in zip(orig_table['FIBERID'], orig_table['MJD']):
	fitFile = join( os.environ[env], dir, plate, "spFly-"+plate.zfill(4)+"-"+str(mjd)+"-"+str(fiber).zfill(4)+suffix)
	if os.path.isfile(fitFile):
		#print fitFile
		table_entry_1, headers_1 = get_table_entry_full( hduSPM=fits.open(fitFile)[1] )
		table_entry_2, headers_2 = get_table_entry_full( hduSPM=fits.open(fitFile)[2] )
		table_entry_3, headers_3 = get_table_entry_full( hduSPM=fits.open(fitFile)[3] )
		headers = headers_1 + headers_2 + headers_3
		table_all.append(n.hstack((table_entry_1, table_entry_2, table_entry_3)))
		#print len(table_all[-1])
		fitFileLast = fitFile
	else:
		table_all.append(n.ones(195)*dV)

newDat = n.transpose(table_all)

all_cols = []
for data_array, head in zip(newDat, headers.split()):
	all_cols.append(fits.Column(name=head, format='D', array=data_array))

new_cols = fits.ColDefs(all_cols)
tbhdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)


sp_cha=fits.open(fitFileLast)[0]

prihdr['model']  = sp_cha.header['model']
prihdr['ageMin'] = sp_cha.header['ageMin']
prihdr['ageMax'] = sp_cha.header['ageMax']
prihdr['Zmin']   = sp_cha.header['Zmin']
prihdr['Zmax']   = sp_cha.header['Zmax']

prihdu = fits.PrimaryHDU(header=prihdr)

hdu = fits.HDUList([prihdu, tbhdu])

if os.path.isfile(path_2_out_file):
    os.remove(path_2_out_file)

hdu.writeto(path_2_out_file)

##################################################
##################################################


test = fits.open(path_2_out_file)

converged3 = (test[1].data['Salpeter_age_lightW'] != dV ) & (test[1].data['Chabrier_age_lightW'] != dV ) & (test[1].data['Kroupa_age_lightW']   != dV )

# now creates the figure per model 
fig = p.figure(0, figsize = (8, 8), frameon=False)#, tight_layout=True)
rect = 0.2, 0.15, 0.85, 0.95
#ax = fig.add_axes(rect, frameon=False)

# panel age distribution
fig.add_subplot(2,2,1)

bins = n.arange(6,10,0.1)
nn_s, bb = n.histogram(n.log10(test[1].data['Salpeter_age_lightW'][converged3]), normed = True, bins=bins)
nn_k, bb = n.histogram(n.log10(test[1].data['Kroupa_age_lightW'][converged3])  , normed = True, bins=bins)
nn_c, bb = n.histogram(n.log10(test[1].data['Chabrier_age_lightW'][converged3]), normed = True, bins=bins)
xb = (bins[:-1]+bins[1:])/2.

p.plot(xb, nn_s, label="Salpeter"  , rasterized =True )
p.plot(xb, nn_k, label="Kroupa"    , rasterized =True )
p.plot(xb, nn_c, label="Chabrier"  , rasterized =True )
p.legend(frameon=False)
p.xlabel('log(age/[yr])')
p.ylabel('Normed cumulative distribution')
p.title("plate=" + plate)
p.grid()

# panel stellar mass 
fig.add_subplot(2,2,2)

bins = n.arange(8,12.5,0.1)
nn_s, bb = n.histogram(test[1].data['Salpeter_stellar_mass'][converged3], normed = True, bins=bins)
nn_k, bb = n.histogram(test[1].data[  'Kroupa_stellar_mass'][converged3], normed = True, bins=bins)
nn_c, bb = n.histogram(test[1].data['Chabrier_stellar_mass'][converged3], normed = True, bins=bins)
xb = (bins[:-1]+bins[1:])/2.

p.plot(xb, nn_s, label="Salpeter"  , rasterized =True )
p.plot(xb, nn_k, label="Kroupa"    , rasterized =True )
p.plot(xb, nn_c, label="Chabrier"  , rasterized =True )

p.xlabel(r'$\log(mass/[M_\odot])$')
#p.ylabel('Normed distribution')
p.grid()

# panels stellar mass difference panels
fig.add_subplot(2,2,3)
p.plot(test[1].data[  'Kroupa_stellar_mass'][converged3], test[1].data['Salpeter_stellar_mass'][converged3] - test[1].data[  'Kroupa_stellar_mass'][converged3], 'k+', rasterized = True, label='Salpeter - Kroupa')
p.xlabel(r'$\log(mass/[M_\odot])$ Kroupa')
p.ylabel(r'$\Delta(\log(M)$')
p.grid()
p.ylim((-0.6,0.6))
p.legend(frameon=False)

fig.add_subplot(2,2,4)
p.plot(test[1].data[  'Kroupa_stellar_mass'][converged3], test[1].data['Chabrier_stellar_mass'][converged3] - test[1].data[  'Kroupa_stellar_mass'][converged3], 'k+', rasterized = True, label='Chabrier - Kroupa')
p.xlabel(r'$\log(mass/[M_\odot])$ Kroupa')
#p.ylabel('mass (Chabrier - Kroupa)')
p.ylim((-0.6,0.6))
p.legend(frameon=False)
p.grid()
p.savefig(path_2_im_file)
p.clf()

print time.time()-t0t


sys.exit()

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

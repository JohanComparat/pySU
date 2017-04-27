import time
t0t=time.time()
from os.path import join
import os
import numpy as n
import glob 
import sys 
import astropy.io.fits as fits

plate   = sys.argv[1]
mjd     = sys.argv[2] 
fiberid = sys.argv[3] 

env = 'EBOSSDR14_DIR'

#dirs = ['stellarpop-m11-salpeter', 'stellarpop-m11-kroupa', 'stellarpop-m11-chabrier', 'stellarpop-m11-salpeter-stelib', 'stellarpop-m11-kroupa-stelib', 'stellarpop-m11-chabrier-stelib', 'stellarpop-m11-salpeter-elodie', 'stellarpop-m11-kroupa-elodie', 'stellarpop-m11-chabrier-elodie'] 
#suffixs = ["-ss.fits", "-kr.fits", "-cha.fits", "-ss.fits", "-kr.fits", "-cha.fits", "-ss.fits", "-kr.fits", "-cha.fits"]

dirs = ['stellarpop-m11-salpeter', 'stellarpop-m11-kroupa', 'stellarpop-m11-chabrier'] 
suffixs = ["-ss.fits", "-kr.fits", "-cha.fits"]

print plate, mjd, fiberid

sp_cha = fits.open(os.path.join(os.environ[env], dirs[0], 'stellarpop', plate, 'spFly-'+plate+'-'+mjd+'-'+fiberid+suffixs[0]))
sp_kr  = fits.open(os.path.join(os.environ[env], dirs[1], 'stellarpop', plate, 'spFly-'+plate+'-'+mjd+'-'+fiberid+suffixs[1]))
sp_sa  = fits.open(os.path.join(os.environ[env], dirs[2], 'stellarpop', plate, 'spFly-'+plate+'-'+mjd+'-'+fiberid+suffixs[2]))
sp_cha_nd = fits.open(os.path.join(os.environ[env],dirs[3], 'stellarpop', plate, 'spFly-'+plate+'-'+mjd+'-'+fiberid+suffixs[3]))
sp_kr_nd = fits.open(os.path.join(os.environ[env], dirs[4], 'stellarpop', plate, 'spFly-'+plate+'-'+mjd+'-'+fiberid+suffixs[4]))
sp_sa_nd = fits.open(os.path.join(os.environ[env], dirs[5], 'stellarpop', plate, 'spFly-'+plate+'-'+mjd+'-'+fiberid+suffixs[5]))
sp_cha_el = fits.open(os.path.join(os.environ[env],dirs[6], 'stellarpop', plate, 'spFly-'+plate+'-'+mjd+'-'+fiberid+suffixs[6]))
sp_kr_el = fits.open(os.path.join(os.environ[env], dirs[7], 'stellarpop', plate, 'spFly-'+plate+'-'+mjd+'-'+fiberid+suffixs[7]))
sp_sa_el = fits.open(os.path.join(os.environ[env], dirs[8], 'stellarpop', plate, 'spFly-'+plate+'-'+mjd+'-'+fiberid+suffixs[8]))

out_dir = os.path.join(os.environ[env], 'stellarpop', plate)
if os.path.isdir(out_dir)==False:
	os.makedirs(out_dir)
out_file = 'spFly-'+plate+'-'+mjd+'-'+fiberid+'.fits'
path_2_out_file = os.path.join(out_dir, out_file)


def create_tbhdu(sp_cha, imf, lib):
    c1 = fits.Column(name='wavelength', format='D', unit='Angstrom', array=sp_cha[1].data['wavelength'])
    c2 = fits.Column(name='model_flux', format='D', unit='1e-17 erg/cm2/s', array=sp_cha[1].data['firefly_model'])
    
    coldefs = fits.ColDefs([c1, c2])
    tbhdu = fits.BinTableHDU.from_columns(coldefs)
    
    tbhdu.header['HIERARCH library'] = lib 
    tbhdu.header['HIERARCH IMF'] = imf 
    tbhdu.header['HIERARCH age_lightW']             = sp_cha[1].header['HIERARCH age_lightW_mean']                                
    tbhdu.header['HIERARCH age_lightW_up']          = sp_cha[1].header['HIERARCH age_lightW_mean_up']                                        
    tbhdu.header['HIERARCH age_lightW_low']         = sp_cha[1].header['HIERARCH age_lightW_mean_low']                                       
    tbhdu.header['HIERARCH metallicity_lightW']     = sp_cha[1].header['HIERARCH metallicity_lightW_mean']                                     
    tbhdu.header['HIERARCH metallicity_lightW_up']  = sp_cha[1].header['HIERARCH metallicity_lightW_mean_up']                                          
    tbhdu.header['HIERARCH metallicity_lightW_low'] = sp_cha[1].header['HIERARCH metallicity_lightW_mean_low']                                           
    tbhdu.header['HIERARCH age_massW']              = sp_cha[1].header['HIERARCH age_massW_mean']                                                 
    tbhdu.header['HIERARCH age_massW_up']           = sp_cha[1].header['HIERARCH age_massW_mean_up']                                                    
    tbhdu.header['HIERARCH age_massW_low']          = sp_cha[1].header['HIERARCH age_massW_mean_low']                                                   
    tbhdu.header['HIERARCH metallicity_massW']      = sp_cha[1].header['HIERARCH metallicity_massW_mean']                                       
    tbhdu.header['HIERARCH metallicity_massW_up']   = sp_cha[1].header['HIERARCH metallicity_massW_mean_up']                                           
    tbhdu.header['HIERARCH metallicity_massW_low']  = sp_cha[1].header['HIERARCH metallicity_massW_mean_low']                                          
    tbhdu.header['HIERARCH EBV']                    = sp_cha[1].header['HIERARCH EBV']                                                                  
    tbhdu.header['HIERARCH stellar_mass']           = sp_cha[1].header['HIERARCH stellar_mass_mean']                                                
    tbhdu.header['HIERARCH stellar_mass_up']        = sp_cha[1].header['HIERARCH stellar_mass_mean_up']                                                   
    tbhdu.header['HIERARCH stellar_mass_low']       = sp_cha[1].header['HIERARCH stellar_mass_mean_low']                                                  
    tbhdu.header['HIERARCH ssp_number']             = sp_cha[1].header['HIERARCH ssp_number']
    
    for el in sp_cha[1].header[33:]:
        tbhdu.header['HIERARCH '+el] = sp_cha[1].header[el]
    
    return tbhdu

tbhdu_cha = create_tbhdu(sp_cha, 'Chabrier', 'MILES')
tbhdu_kr  = create_tbhdu(sp_kr, 'Kroupa'   , 'MILES')
tbhdu_sa  = create_tbhdu(sp_sa, 'Salpeter' , 'MILES')
tbhdu_cha_nd = create_tbhdu(sp_cha_nd, 'Chabrier', 'STELIB')
tbhdu_kr_nd  = create_tbhdu(sp_kr_nd, 'Kroupa'   , 'STELIB')
tbhdu_sa_nd  = create_tbhdu(sp_sa_nd, 'Salpeter' , 'STELIB')
tbhdu_cha_el = create_tbhdu(sp_cha_el, 'Chabrier', 'ELODIE')
tbhdu_kr_el  = create_tbhdu(sp_kr_el, 'Kroupa'   , 'ELODIE')
tbhdu_sa_el  = create_tbhdu(sp_sa_el, 'Salpeter' , 'ELODIE')

prihdr = fits.Header()

prihdr['file']   = out_file
prihdr['plate']  = int(plate)
prihdr['mjd']    = int(mjd)
prihdr['fiberid']= int(fiberid)

prihdr['fitter'] = 'FIREFLY'
prihdr['model']  = sp_kr[0].header['model']
prihdr['ageMin'] = sp_kr[0].header['ageMin']
prihdr['ageMax'] = sp_kr[0].header['ageMax']
prihdr['Zmin']   = sp_kr[0].header['Zmin']
prihdr['Zmax']   = sp_kr[0].header['Zmax']

prihdr['HIERARCH age_universe'] = sp_cha[1].header['HIERARCH age_universe']                                             
prihdr['HIERARCH redshift']     = sp_cha[1].header['HIERARCH redshift']                                         

prihdu = fits.PrimaryHDU(header=prihdr)

thdulist = fits.HDUList([prihdu, tbhdu_cha, tbhdu_kr, tbhdu_sa, tbhdu_cha_nd, tbhdu_kr_nd, tbhdu_sa_nd, tbhdu_cha_el, tbhdu_kr_el, tbhdu_sa_el ])
if os.path.isfile(path_2_out_file ):
	os.remove(path_2_out_file )

thdulist.writeto( path_2_out_file )

print time.time()-t0t


sys.exit()









init_cat = join( os.environ['EBOSSDR14_DIR'], "catalogs", "perPlate", "sp-"+plate.zfill(4)+".fits")
plate_catalog = join( os.environ['EBOSSDR14_DIR'], dir, "catalogs", "spFly-"+plate.zfill(4)+".fits")

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

#aaa=orig_cols.del_col('CHUNK'                      )
#aaa=orig_cols.del_col('PLATEQUALITY'           )
#aaa=orig_cols.del_col('PLATESN2'                  )
aaa=orig_cols.del_col('DEREDSN2'                )
aaa=orig_cols.del_col('LAMBDA_EFF'             )
aaa=orig_cols.del_col('BLUEFIBER'                )
aaa=orig_cols.del_col('ZOFFSET'                   )
aaa=orig_cols.del_col('SPECPRIMARY'            )
aaa=orig_cols.del_col('SPECBOSS'                 )
#aaa=orig_cols.del_col('BOSS_SPECOBJ_ID'    )

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

#aaa=orig_cols.del_col('SN_MEDIAN_ALL'       )
#aaa=orig_cols.del_col('SN_MEDIAN'               )
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


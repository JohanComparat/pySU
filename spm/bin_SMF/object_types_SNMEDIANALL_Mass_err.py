"""
Plots log stellar mass vs. log(1+z) for each PROGRAMNAME.

"""

import os
import numpy as n
import astropy.io.fits as fits
import matplotlib
matplotlib.rcParams['agg.path.chunksize'] = 20000000
matplotlib.rcParams.update({'font.size': 13})
#matplotlib.use('Agg')
import matplotlib.pyplot as p


out_dir = os.path.join(os.environ['HOME'], 'wwwDir/firefly_catalogs/figures')

imfs = ["Chabrier_ELODIE_", "Chabrier_MILES_", "Chabrier_STELIB_", "Kroupa_ELODIE_", "Kroupa_MILES_", "Kroupa_STELIB_",  "Salpeter_ELODIE_", "Salpeter_MILES_", "Salpeter_STELIB_" ]

sdss_dir = os.path.join(os.environ['OBS_REPO'], 'SDSS', 'dr14')
path_2_sdss_cat = os.path.join( sdss_dir, 'firefly', "FireflyGalaxySdss26.fits" )
path_2_eboss_cat = os.path.join( sdss_dir, 'firefly', "FireflyGalaxyEbossDR14.fits" )

# DEEP SURVEYS
deep2_dir = os.path.join(os.environ['OBS_REPO'], 'DEEP2')
path_2_deep2_cat = os.path.join( deep2_dir, "zcat.deep2.dr4.v4.LFcatalogTC.Planck15.spm.v3.fits" )


# path_2_F16_cat = os.path.join( sdss_dir, "RA_DEC_z_w_fluxOII_Mstar_grcol_Mr_lumOII.dat" )

# OPENS THE CATALOGS
deep2   = fits.open(path_2_deep2_cat)[1].data
sdss   = fits.open(path_2_sdss_cat)[1].data
boss   = fits.open(path_2_eboss_cat)[1].data

redshift_reliable_boss =  (boss['CLASS_NOQSO'] == "GALAXY") & ( boss['Z_ERR_NOQSO'] > 0.0) & (boss['ZWARNING_NOQSO'] == 0) & (boss['Z_NOQSO']>0.001) & (boss['Z_NOQSO'] > boss['Z_ERR_NOQSO'] ) # (boss['SN_MEDIAN_ALL'] > 0.1 ) & 
redshift_reliable_sdss =  (sdss['CLASS'] == "GALAXY")       & ( sdss['Z_ERR'] > 0.0)       & (sdss['ZWARNING'] == 0)       & (sdss['Z'] > 0.001) & (sdss['Z'] > sdss['Z_ERR'] ) # (sdss['SN_MEDIAN_ALL'] > 0.1 ) &


imf = imfs[0]
stellar_mass = imf+'stellar_mass'

error_reliable_boss = (boss[stellar_mass+'_up'] > boss[stellar_mass+'_low'] ) & (boss[stellar_mass+'_up'] > 0. ) & ( boss[stellar_mass+'_low'] > 1e6 ) & (boss[stellar_mass+'_up'] < 1e14 ) & ( boss[stellar_mass+'_low'] < 1e14 ) 
error_reliable_sdss = (sdss[stellar_mass+'_up'] > sdss[stellar_mass+'_low'] ) & (sdss[stellar_mass+'_up'] > 0. ) & ( sdss[stellar_mass+'_low'] > 1e6 ) & (sdss[stellar_mass+'_up'] < 1e14 ) & ( sdss[stellar_mass+'_low'] < 1e14 ) 

mass_reliable_boss_02 = (boss[stellar_mass] > 1e6 ) & ( boss[stellar_mass] < 1e14 ) & ((n.log10(boss[stellar_mass+'_up']) - n.log10(boss[stellar_mass+'_low']))/2. < 0.2 )
mass_reliable_sdss_02 = (sdss[stellar_mass] > 1e6 ) & ( sdss[stellar_mass] < 1e14 ) & ((n.log10(sdss[stellar_mass+'_up']) - n.log10(sdss[stellar_mass+'_low']))/2. < 0.2 )
mass_reliable_boss_04 = (boss[stellar_mass] > 1e6 ) & ( boss[stellar_mass] < 1e14 ) & ((n.log10(boss[stellar_mass+'_up']) - n.log10(boss[stellar_mass+'_low']))/2. < 0.4 )
mass_reliable_sdss_04 = (sdss[stellar_mass] > 1e6 ) & ( sdss[stellar_mass] < 1e14 ) & ((n.log10(sdss[stellar_mass+'_up']) - n.log10(sdss[stellar_mass+'_low']))/2. < 0.4 )

ok_boss = (error_reliable_boss) & (redshift_reliable_boss) & ( boss['SN_MEDIAN_ALL']>0)# & (mass_reliable_boss_02) 
ok_sdss = (error_reliable_sdss) & (redshift_reliable_sdss) & ( sdss['SN_MEDIAN_ALL']>0)# & (mass_reliable_sdss_02) 

####################SN MEDIAN ALL vs. MASS ERROR ##################

zz = n.hstack(( boss['Z_NOQSO'][ok_boss], sdss['Z'][ok_sdss]))

Ms_log_err = n.hstack(( (n.log10(boss[stellar_mass+'_up'][ok_boss]) - n.log10(boss[stellar_mass+'_low'][ok_boss]))/2., (n.log10(sdss[stellar_mass+'_up'][ok_sdss]) - n.log10(sdss[stellar_mass+'_low'][ok_sdss]))/2.))

sn_all = n.hstack(( boss['SN_MEDIAN_ALL'][ok_boss], sdss['SN_MEDIAN_ALL'][ok_sdss] ))
sn_u = n.hstack(( boss['SN_MEDIAN'].T[0][ok_boss], sdss['SN_MEDIAN'].T[0][ok_sdss] ))
sn_g = n.hstack(( boss['SN_MEDIAN'].T[1][ok_boss], sdss['SN_MEDIAN'].T[1][ok_sdss] ))
sn_r = n.hstack(( boss['SN_MEDIAN'].T[2][ok_boss], sdss['SN_MEDIAN'].T[2][ok_sdss] ))
sn_i = n.hstack(( boss['SN_MEDIAN'].T[3][ok_boss], sdss['SN_MEDIAN'].T[3][ok_sdss] ))
sn_z = n.hstack(( boss['SN_MEDIAN'].T[4][ok_boss], sdss['SN_MEDIAN'].T[4][ok_sdss] ))

def plot_it(z_min, sn_all, ssn="SN_MEDIAN_ALL_"):
	oo = (zz>z_min)&(zz<z_min+0.1)
	p.figure(1, (4.5, 4.5))
	p.axes([0.2,0.2,0.7,0.7])
	p.plot(Ms_log_err[oo], sn_all[oo],  linestyle='None', color='k', marker=',', label='z='+str(n.round(z_min+0.05,2)), rasterized=True, alpha=0.1 )
	err_bins = n.arange(0,0.81,0.05)
	for err_min in err_bins:
		ee = (Ms_log_err[oo]>err_min)&(Ms_log_err[oo]<err_min+0.025)
		mean = n.median(sn_all[oo][ee])
		std = n.std(sn_all[oo][ee])
		print([err_min+0.025, err_min+0.025], [mean+std, mean-std])
		p.plot([err_min+0.025, err_min+0.025], [mean+std, mean-std], 'r', lw=1)
	

	p.axhline(20, ls='dashed', color='b', label='SN=20, 5')
	p.axhline(5, ls='dashed', color='b')
	p.ylabel("SN MEDIAN "+ssn.split('_')[-2])
	p.xlabel(r"$\sigma_M$ [dex]")
	#p.xlabel(r'$\log_{10}(1+z)$')
	p.yscale('log')
	p.legend(loc=0, frameon = False)
	p.xlim((0., 0.6))
	p.ylim((0.5, 60))
	p.grid()
	#p.title('SDSS and BOSS')
	p.savefig(os.path.join(out_dir, ssn+"stellar_mass_err_sdss_boss"+str(int(z_min*10)).zfill(2)+".jpg" ))
	p.clf()

for z_min in n.arange(0, 1.1, 0.1):
	plot_it(z_min, sn_all, ssn="SN_MEDIAN_ALL_")
	plot_it(z_min, sn_u, ssn="SN_MEDIAN_U_")
	plot_it(z_min, sn_g, ssn="SN_MEDIAN_G_")
	plot_it(z_min, sn_r, ssn="SN_MEDIAN_R_")
	plot_it(z_min, sn_i, ssn="SN_MEDIAN_I_")
	plot_it(z_min, sn_z, ssn="SN_MEDIAN_Z_")
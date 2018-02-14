"""
Plots log stellar mass vs. log(1+z) for each PROGRAMNAME.

"""
from lib_spm import *

#import matplotlib
#matplotlib.rcParams['agg.path.chunksize'] = 2000000
#matplotlib.rcParams.update({'font.size': 13})
##matplotlib.use('Agg')
#import matplotlib.pyplot as p

out_dir = os.path.join(os.environ['OBS_REPO'], 'spm', 'results', 'mass-snr')

imfs = ["Chabrier_ELODIE_", "Chabrier_MILES_", "Chabrier_STELIB_", "Kroupa_ELODIE_", "Kroupa_MILES_", "Kroupa_STELIB_",  "Salpeter_ELODIE_", "Salpeter_MILES_", "Salpeter_STELIB_" ]

#sdss_dir = os.path.join(os.environ['OBS_REPO'], 'SDSS', 'dr14')
#path_2_sdss_cat = os.path.join( sdss_dir, 'firefly', "FireflyGalaxySdss26.fits" )
#path_2_eboss_cat = os.path.join( sdss_dir, 'firefly', "FireflyGalaxyEbossDR14.fits" )

## DEEP SURVEYS
#deep2_dir = os.path.join(os.environ['OBS_REPO'], 'DEEP2')
#path_2_deep2_cat = os.path.join( deep2_dir, "zcat.deep2.dr4.v4.LFcatalogTC.Planck15.spm.v3.fits" )


## path_2_F16_cat = os.path.join( sdss_dir, "RA_DEC_z_w_fluxOII_Mstar_grcol_Mr_lumOII.dat" )

## OPENS THE CATALOGS
#deep2   = fits.open(path_2_deep2_cat)[1].data
#sdss   = fits.open(path_2_sdss_cat)[1].data
#boss   = fits.open(path_2_eboss_cat)[1].data

redshift_reliable_boss =  (boss['CLASS_NOQSO'] == "GALAXY") & ( boss['Z_ERR_NOQSO'] > 0.0) & (boss['ZWARNING_NOQSO'] == 0) & (boss['Z_NOQSO']>0.001) & (boss['Z_NOQSO'] > boss['Z_ERR_NOQSO'] ) # (boss['SNR_ALL'] > 0.1 ) & 
redshift_reliable_sdss =  (sdss['CLASS'] == "GALAXY")       & ( sdss['Z_ERR'] > 0.0)       & (sdss['ZWARNING'] == 0)       & (sdss['Z'] > 0.001) & (sdss['Z'] > sdss['Z_ERR'] ) # (sdss['SNR_ALL'] > 0.1 ) &


imf = imfs[0]
stellar_mass = imf+'stellar_mass'

error_reliable_boss = (boss[stellar_mass+'_up_1sig'] > boss[stellar_mass+'_low_1sig'] ) & (boss[stellar_mass+'_up_1sig'] > 0. ) & ( boss[stellar_mass+'_low_1sig'] > 0. ) & (boss[stellar_mass+'_up_1sig'] < 1e14 ) & ( boss[stellar_mass+'_low_1sig'] < 1e14 ) 
error_reliable_sdss = (sdss[stellar_mass+'_up_1sig'] > sdss[stellar_mass+'_low_1sig'] ) & (sdss[stellar_mass+'_up_1sig'] > 0. ) & ( sdss[stellar_mass+'_low_1sig'] > 0. ) & (sdss[stellar_mass+'_up_1sig'] < 1e14 ) & ( sdss[stellar_mass+'_low_1sig'] < 1e14 ) 

mass_reliable_boss_02 = (boss[stellar_mass] > 1e6 ) & ( boss[stellar_mass] < 1e14 ) & ((n.log10(boss[stellar_mass+'_up_1sig']) - n.log10(boss[stellar_mass+'_low_1sig']))/2. < 0.2 )
mass_reliable_sdss_02 = (sdss[stellar_mass] > 1e6 ) & ( sdss[stellar_mass] < 1e14 ) & ((n.log10(sdss[stellar_mass+'_up_1sig']) - n.log10(sdss[stellar_mass+'_low_1sig']))/2. < 0.2 )
mass_reliable_boss_04 = (boss[stellar_mass] > 1e6 ) & ( boss[stellar_mass] < 1e14 ) & ((n.log10(boss[stellar_mass+'_up_1sig']) - n.log10(boss[stellar_mass+'_low_1sig']))/2. < 0.4 )
mass_reliable_sdss_04 = (sdss[stellar_mass] > 1e6 ) & ( sdss[stellar_mass] < 1e14 ) & ((n.log10(sdss[stellar_mass+'_up_1sig']) - n.log10(sdss[stellar_mass+'_low_1sig']))/2. < 0.4 )

ok_boss_02 = (error_reliable_boss) & (mass_reliable_boss_02) & (redshift_reliable_boss) & ( boss['SNR_ALL']>0)
ok_sdss_02 = (error_reliable_sdss) & (mass_reliable_sdss_02) & (redshift_reliable_sdss) & ( sdss['SNR_ALL']>0)
ok_boss_04 = (error_reliable_boss) & (mass_reliable_boss_04) & (redshift_reliable_boss) & ( boss['SNR_ALL']>0)
ok_sdss_04 = (error_reliable_sdss) & (mass_reliable_sdss_04) & (redshift_reliable_sdss) & ( sdss['SNR_ALL']>0)
print( "boss 02",len(ok_boss_02.nonzero()[0]))
print( "sdss 02",len(ok_sdss_02.nonzero()[0]))
print( "boss 04",len(ok_boss_04.nonzero()[0]))
print( "sdss 04",len(ok_sdss_04.nonzero()[0]))

zz_02 = n.hstack(( boss['Z_NOQSO'][ok_boss_02], sdss['Z'][ok_sdss_02]))
sn_02 = n.hstack(( boss['SNR_ALL'][ok_boss_02], sdss['SNR_ALL'][ok_sdss_02] ))

zz_04 = n.hstack(( boss['Z_NOQSO'][ok_boss_04], sdss['Z'][ok_sdss_04]))
sn_04 = n.hstack(( boss['SNR_ALL'][ok_boss_04], sdss['SNR_ALL'][ok_sdss_04] ))


####################SN MEDIAN ALL vs. MASS ERROR ##################


Ms_log_err = n.hstack(( (n.log10(boss[stellar_mass+'_up_1sig'][redshift_reliable_boss & error_reliable_boss]) - n.log10(boss[stellar_mass+'_low_1sig'][redshift_reliable_boss & error_reliable_boss]))/2., (n.log10(sdss[stellar_mass+'_up_1sig'][redshift_reliable_sdss & error_reliable_sdss]) - n.log10(sdss[stellar_mass+'_low_1sig'][redshift_reliable_sdss & error_reliable_sdss]))/2.))

sn_all = n.hstack(( boss['SNR_ALL'][redshift_reliable_boss & error_reliable_boss], sdss['SNR_ALL'][redshift_reliable_sdss & error_reliable_sdss] ))

p.figure(1, (4.5, 4.5))
p.axes([0.2,0.2,0.7,0.7])
#p.plot(zz_04, Ms_04, 'r,', rasterized=True, alpha=0.5, label=r' $\sigma_M < 0.4$')
p.plot(Ms_log_err, sn_all, 'k,', rasterized=True, alpha=0.5)
err_bins = n.arange(0,1.5,0.05)
for err_min in err_bins:
	oo = (Ms_log_err>err_min)&(Ms_log_err<err_min+0.025)
	mean = n.median(sn_all[oo])
	std = n.std(sn_all[oo])
	print([err_min+0.05, err_min+0.05], [mean+std, mean-std])
	p.plot([err_min+0.05, err_min+0.05], [mean+std, mean-std], 'r', lw=1)
	
p.ylabel("SN MEDIAN ALL")
#p.axvline(n.log10(3.), ls='dashed', label='z=0.1, 0.5, 1, 2')
#p.axvline(n.log10(2.), ls='dashed')#, label='z=1')
#p.axvline(n.log10(1.1),ls='dashed')#, label='z=0.1')
#p.axvline(n.log10(1.5),ls='dashed')#, label='z=0.5')
p.xlabel(r"$\sigma_M$ [dex]")
#p.xlabel(r'$\log_{10}(1+z)$')
p.yscale('log')
p.legend(loc=0, frameon = False)
p.xlim((0., 2.))
p.ylim((0.1, 200))
p.grid()
p.title('SDSS and BOSS')
p.savefig(os.path.join(out_dir, "SNR_ALL_stellar_mass_err_sdss_boss.png" ))
p.clf()

#################### Z vs. SN MEDIAN ALL ##########################
p.figure(1, (4.5, 4.5))
p.axes([0.2,0.2,0.7,0.7])
#p.plot(zz_04, Ms_04, 'r,', rasterized=True, alpha=0.5, label=r' $\sigma_M < 0.4$')
p.plot(zz_02, sn_02, 'k,', rasterized=True, alpha=0.5, label=r' $\sigma_M < 0.2$')
z_bins = n.arange(0,1.4,0.1)
for z_min in z_bins:
	oo = (zz_02>z_min)&(zz_02<z_min+0.1)
	mean = n.median(sn_02[oo])
	std = n.std(sn_02[oo])
	print([z_min+0.05, z_min+0.05], [mean+std, mean-std])
	p.plot([z_min+0.05, z_min+0.05], [mean+std, mean-std], 'r', lw=1)
	
p.ylabel("SN MEDIAN ALL")
#p.axvline(n.log10(3.), ls='dashed', label='z=0.1, 0.5, 1, 2')
#p.axvline(n.log10(2.), ls='dashed')#, label='z=1')
#p.axvline(n.log10(1.1),ls='dashed')#, label='z=0.1')
#p.axvline(n.log10(1.5),ls='dashed')#, label='z=0.5')
p.xlabel('redshift')
#p.xlabel(r'$\log_{10}(1+z)$')
p.yscale('log')
p.legend(loc=0, frameon = False)
p.xlim((0.0, 1.4))
p.ylim((0.1, 200))
p.grid()
p.title('SDSS and BOSS')
p.savefig(os.path.join(out_dir, "SNR_ALL_redshift_sdss_boss.png" ))
p.clf()



p.figure(1, (4.5, 4.5))
p.axes([0.2,0.2,0.7,0.7])
p.plot(zz_04, sn_04, 'k,', rasterized=True, alpha=0.5, label=r' $\sigma_M < 0.4$')
#p.plot(zz_02, Ms_02, 'k,', rasterized=True, alpha=0.5, label=r' $\sigma_M < 0.2$')
for z_min in z_bins:
	oo = (zz_04>z_min)&(zz_04<z_min+0.1)
	mean = n.median(sn_04[oo])
	std = n.std(sn_04[oo])
	print([z_min+0.05, z_min+0.05], [mean+std, mean-std])
	p.plot([z_min+0.05, z_min+0.05], [mean+std, mean-std], 'r', lw=1)
	
p.ylabel("SN MEDIAN ALL")
#p.axvline(n.log10(3.), ls='dashed', label='z=0.1, 0.5, 1, 2')
#p.axvline(n.log10(2.), ls='dashed')#, label='z=1')
#p.axvline(n.log10(1.1),ls='dashed')#, label='z=0.1')
#p.axvline(n.log10(1.5),ls='dashed')#, label='z=0.5')
p.xlabel('redshift')
#p.xlabel(r'$\log_{10}(1+z)$')
p.yscale('log')
p.legend(loc=0, frameon = False)
p.xlim((0.0, 1.4))
p.ylim((0.1, 200))
p.grid()
p.title('SDSS and BOSS')
p.savefig(os.path.join(out_dir, "SNR_ALL_redshift_sdss_boss_04.png" ))
p.clf()



p.figure(1, (4.5, 4.5))
p.axes([0.2,0.2,0.7,0.7])
#p.plot(zz_04, Ms_04, 'r,', rasterized=True, alpha=0.5, label=r' $\sigma_M < 0.4$')
p.plot(zz_02, Ms_02, 'k,', rasterized=True, alpha=0.5, label=r' $\sigma_M < 0.2$')
z_bins = n.arange(0,1.4,0.1)
for z_min in z_bins:
	oo = (zz_02>z_min)&(zz_02<z_min+0.1)
	mean = n.median(Ms_02[oo])
	std = n.std(Ms_02[oo])
	print([z_min+0.05, z_min+0.05], [mean+std, mean-std])
	p.plot([z_min+0.05, z_min+0.05], [mean+std, mean-std], 'r', lw=1)
	
p.ylabel("SN MEDIAN ALL")
#p.axvline(n.log10(3.), ls='dashed', label='z=0.1, 0.5, 1, 2')
#p.axvline(n.log10(2.), ls='dashed')#, label='z=1')
#p.axvline(n.log10(1.1),ls='dashed')#, label='z=0.1')
#p.axvline(n.log10(1.5),ls='dashed')#, label='z=0.5')
p.xlabel('redshift')
#p.xlabel(r'$\log_{10}(1+z)$')
p.yscale('log')
p.legend(loc=0, frameon = False)
p.xlim((0.0, 1.4))
p.ylim((0.1, 200))
p.grid()
p.title('SDSS and BOSS')
p.savefig(os.path.join(out_dir, "SNR_ALL_redshift_sdss_boss.png" ))
p.clf()



p.figure(1, (4.5, 4.5))
p.axes([0.2,0.2,0.7,0.7])
p.plot(zz_04, Ms_04, 'k,', rasterized=True, alpha=0.5, label=r' $\sigma_M < 0.4$')
#p.plot(zz_02, Ms_02, 'k,', rasterized=True, alpha=0.5, label=r' $\sigma_M < 0.2$')
for z_min in z_bins:
	oo = (zz_04>z_min)&(zz_04<z_min+0.1)
	mean = n.median(Ms_04[oo])
	std = n.std(Ms_04[oo])
	print([z_min+0.05, z_min+0.05], [mean+std, mean-std])
	p.plot([z_min+0.05, z_min+0.05], [mean+std, mean-std], 'r', lw=1)
	
p.ylabel("SN MEDIAN ALL")
#p.axvline(n.log10(3.), ls='dashed', label='z=0.1, 0.5, 1, 2')
#p.axvline(n.log10(2.), ls='dashed')#, label='z=1')
#p.axvline(n.log10(1.1),ls='dashed')#, label='z=0.1')
#p.axvline(n.log10(1.5),ls='dashed')#, label='z=0.5')
p.xlabel('redshift')
#p.xlabel(r'$\log_{10}(1+z)$')
p.yscale('log')
p.legend(loc=0, frameon = False)
p.xlim((0.0, 1.4))
p.ylim((0.1, 200))
p.grid()
p.title('SDSS and BOSS')
p.savefig(os.path.join(out_dir, "SNR_ALL_redshift_sdss_boss_04.png" ))
p.clf()


z_flg = 'ZQUALITY'
z_name = 'ZBEST'
deep2_zOk = (deep2[z_name] > 0.001) & (deep2[z_flg]>=2.) & (deep2[z_name] < 1.7) & (deep2['SSR']>0) & (deep2['TSR']>0) & (deep2['SSR']<=1.0001) & (deep2['TSR']<=1.0001)
deep2_sel_02 = (deep2_zOk) & (deep2[stellar_mass] < 10**14. ) & (deep2[stellar_mass] > 0. )  & (deep2[stellar_mass] >= deep2[stellar_mass+'_low_1sig'] ) & (deep2[stellar_mass] <= deep2[stellar_mass+'_up_1sig'] ) & ( - n.log10(deep2[stellar_mass+'_low_1sig'])  + n.log10(deep2[stellar_mass+'_up_1sig']) < 0.4 )
deep2_sel_04 = (deep2_zOk) & (deep2[stellar_mass] < 10**14. ) & (deep2[stellar_mass] > 0. )  & (deep2[stellar_mass] >= deep2[stellar_mass+'_low_1sig'] ) & (deep2[stellar_mass] <= deep2[stellar_mass+'_up_1sig'] ) & ( - n.log10(deep2[stellar_mass+'_low_1sig'])  + n.log10(deep2[stellar_mass+'_up_1sig']) < 0.8 )

Ms_02 = deep2['SNR_ALL'][deep2_sel_02]
zz_02 = deep2[z_name][deep2_sel_02]
Ms_04 = deep2['SNR_ALL'][deep2_sel_04]
zz_04 = deep2[z_name][deep2_sel_04]

w_deep2 = 1. / (deep2['TSR'] * deep2['SSR'])

print( "deep2 02",len(deep2_sel_02.nonzero()[0]))
print( "deep2 04",len(deep2_sel_04.nonzero()[0]))

p.figure(1, (4.5, 4.5))
p.axes([0.2,0.2,0.7,0.7])
#p.plot(zz_04, Ms_04, 'r+', rasterized=True, label=r' $\sigma_M < 0.4$', alpha=0.5)
p.plot(zz_02, Ms_02, 'kx', rasterized=True, label=r' $\sigma_M < 0.2$', alpha=0.5)
p.ylabel("SN MEDIAN ALL")
#p.axvline(n.log10(3.), ls='dashed', label='z=0.1, 0.5, 1, 2')
#p.axvline(n.log10(2.), ls='dashed')#, label='z=1')
#p.axvline(n.log10(1.1),ls='dashed')#, label='z=0.1')
#p.axvline(n.log10(1.5),ls='dashed')#, label='z=0.5')
p.xlabel('redshift')
#p.xlabel(r'$\log_{10}(1+z)$')
p.yscale('log')
p.legend(loc=0, frameon = False)
p.xlim((0.0, 1.4))
p.ylim((0.1, 200))
p.grid()
p.title('DEEP2')
p.savefig(os.path.join(out_dir, "SNR_ALL_mass_deep2.png" ))
p.clf()

p.figure(1, (4.5, 4.5))
p.axes([0.2,0.2,0.7,0.7])
p.plot(zz_04, Ms_04, 'k+', rasterized=True, label=r' $\sigma_M < 0.4$', alpha=0.5)
#p.plot(zz_02, Ms_02, 'kx', rasterized=True, label=r' $\sigma_M < 0.2$', alpha=0.5)
p.ylabel("SN MEDIAN ALL")
#p.axvline(n.log10(3.), ls='dashed', label='z=0.1, 0.5, 1, 2')
#p.axvline(n.log10(2.), ls='dashed')#, label='z=1')
#p.axvline(n.log10(1.1),ls='dashed')#, label='z=0.1')
#p.axvline(n.log10(1.5),ls='dashed')#, label='z=0.5')
p.xlabel('redshift')
#p.xlabel(r'$\log_{10}(1+z)$')
p.yscale('log')
p.legend(loc=0, frameon = False)
p.xlim((0.0, 1.4))
p.ylim((0.1, 200))
p.grid()
p.title('DEEP2')
p.savefig(os.path.join(out_dir, "SNR_ALL_mass_deep2_04.png" ))
p.clf()




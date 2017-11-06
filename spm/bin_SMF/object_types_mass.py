"""
Plots log stellar mass vs. log(1+z) for each PROGRAMNAME.

"""

import os
import numpy as n
import astropy.io.fits as fits
import matplotlib
matplotlib.rcParams['agg.path.chunksize'] = 2000000
matplotlib.rcParams.update({'font.size': 13})
#matplotlib.use('Agg')
import matplotlib.pyplot as p


out_dir = os.path.join(os.environ['OBS_REPO'], 'spm', 'results', 'mass-redshift-presentation')

imfs = ["Chabrier_ELODIE_", "Chabrier_MILES_", "Chabrier_STELIB_", "Kroupa_ELODIE_", "Kroupa_MILES_", "Kroupa_STELIB_",  "Salpeter_ELODIE_", "Salpeter_MILES_", "Salpeter_STELIB_" ]

cosmos_dir = os.path.join(os.environ['OBS_REPO'], 'COSMOS', 'catalogs' )
path_2_cosmos_cat = os.path.join( cosmos_dir, "photoz-2.0", "photoz_vers2.0_010312.fits")
#path_2_cosmos_cat = os.path.join( cosmos_dir, "COSMOS2015_Laigle+_v1.1.fits.gz")

sdss_dir = os.path.join(os.environ['OBS_REPO'], 'SDSS', 'dr14')
path_2_sdss_cat = os.path.join( sdss_dir, 'firefly', "FireflyGalaxySdss26.fits" )
path_2_eboss_cat = os.path.join( sdss_dir, 'firefly', "FireflyGalaxyEbossDR14.fits" )

# DEEP SURVEYS
deep2_dir = os.path.join(os.environ['OBS_REPO'], 'DEEP2')
path_2_deep2_cat = os.path.join( deep2_dir, "zcat.deep2.dr4.v4.LFcatalogTC.Planck15.spm.v2.fits" )

vipers_dir = os.path.join(os.environ['OBS_REPO'], 'VIPERS')
path_2_vipers_cat = os.path.join( vipers_dir, "VIPERS_W14_summary_v2.1.linesFitted.spm.fits" )

vvds_dir = os.path.join(os.environ['OBS_REPO'], 'VVDS')
path_2_vvdsW_cat = os.path.join( vvds_dir, "catalogs", "VVDS_WIDE_summary.v1.spm.fits" )
path_2_vvdsD_cat = os.path.join( vvds_dir, "catalogs", "VVDS_DEEP_summary.v1.spm.fits" )

# path_2_F16_cat = os.path.join( sdss_dir, "RA_DEC_z_w_fluxOII_Mstar_grcol_Mr_lumOII.dat" )

# OPENS THE CATALOGS
deep2   = fits.open(path_2_deep2_cat)[1].data
sdss   = fits.open(path_2_sdss_cat)[1].data
boss   = fits.open(path_2_eboss_cat)[1].data
cosmos = fits.open(path_2_cosmos_cat)[1].data


imf = imfs[0]
stellar_mass = imf+'stellar_mass'

redshift_reliable_boss =  (boss['CLASS_NOQSO'] == "GALAXY") & ( boss['Z_ERR_NOQSO'] > 0.0) & (boss['ZWARNING_NOQSO'] == 0) & (boss['Z_NOQSO']>0.001) & (boss['Z_NOQSO'] > boss['Z_ERR_NOQSO'] ) # (boss['SN_MEDIAN_ALL'] > 0.1 ) & 
redshift_reliable_sdss =  (sdss['CLASS'] == "GALAXY")       & ( sdss['Z_ERR'] > 0.0)       & (sdss['ZWARNING'] == 0)       & (sdss['Z'] > 0.001) & (sdss['Z'] > sdss['Z_ERR'] ) # (sdss['SN_MEDIAN_ALL'] > 0.1 ) &

error_reliable_boss = (boss[stellar_mass+'_up'] > boss[stellar_mass+'_low'] ) & (boss[stellar_mass+'_up'] > 0. ) & ( boss[stellar_mass+'_low'] > 0. ) & (boss[stellar_mass+'_up'] < 1e14 ) & ( boss[stellar_mass+'_low'] < 1e14 ) 
error_reliable_sdss = (sdss[stellar_mass+'_up'] > sdss[stellar_mass+'_low'] ) & (sdss[stellar_mass+'_up'] > 0. ) & ( sdss[stellar_mass+'_low'] > 0. ) & (sdss[stellar_mass+'_up'] < 1e14 ) & ( sdss[stellar_mass+'_low'] < 1e14 ) 

mass_reliable_boss_02 = (boss[stellar_mass] > 1e6 ) & ( boss[stellar_mass] < 1e14 ) & ((n.log10(boss[stellar_mass+'_up']) - n.log10(boss[stellar_mass+'_low']))/2. < 0.2 )
mass_reliable_sdss_02 = (sdss[stellar_mass] > 1e6 ) & ( sdss[stellar_mass] < 1e14 ) & ((n.log10(sdss[stellar_mass+'_up']) - n.log10(sdss[stellar_mass+'_low']))/2. < 0.2 )
mass_reliable_boss_04 = (boss[stellar_mass] > 1e6 ) & ( boss[stellar_mass] < 1e14 ) & ((n.log10(boss[stellar_mass+'_up']) - n.log10(boss[stellar_mass+'_low']))/2. < 0.4 )
mass_reliable_sdss_04 = (sdss[stellar_mass] > 1e6 ) & ( sdss[stellar_mass] < 1e14 ) & ((n.log10(sdss[stellar_mass+'_up']) - n.log10(sdss[stellar_mass+'_low']))/2. < 0.4 )

ok_boss_02 = (error_reliable_boss) & (mass_reliable_boss_02) & (redshift_reliable_boss)
ok_sdss_02 = (error_reliable_sdss) & (mass_reliable_sdss_02) & (redshift_reliable_sdss)
ok_boss_04 = (error_reliable_boss) & (mass_reliable_boss_04) & (redshift_reliable_boss)
ok_sdss_04 = (error_reliable_sdss) & (mass_reliable_sdss_04) & (redshift_reliable_sdss)
print( "boss 02",len(ok_boss_02.nonzero()[0]))
print( "sdss 02",len(ok_sdss_02.nonzero()[0]))
print( "boss 04",len(ok_boss_04.nonzero()[0]))
print( "sdss 04",len(ok_sdss_04.nonzero()[0]))

zz_02 = n.hstack(( boss['Z_NOQSO'][ok_boss_02], sdss['Z'][ok_sdss_02]))
Ms_02 = n.hstack(( n.log10(boss[stellar_mass][ok_boss_02]), n.log10(sdss[stellar_mass][ok_sdss_02]) ))
zz_02_boss = boss['Z_NOQSO'][ok_boss_02]
zz_02_sdss = sdss['Z'][ok_sdss_02]

zz_04 = n.hstack(( boss['Z_NOQSO'][ok_boss_04], sdss['Z'][ok_sdss_04]))
Ms_04 = n.hstack(( n.log10(boss[stellar_mass][ok_boss_04]), n.log10(sdss[stellar_mass][ok_sdss_04]) ))
zz_04_boss = boss['Z_NOQSO'][ok_boss_04]
zz_04_sdss = sdss['Z'][ok_sdss_04]

p.figure(1, (4.5, 4.5))
p.axes([0.2,0.2,0.7,0.7])
#p.plot(zz_04, Ms_04, 'r,', rasterized=True, alpha=0.5, label=r' $\sigma_M < 0.4$')
p.plot(zz_02, Ms_02, 'k,', rasterized=True, alpha=0.5, label=r' $\sigma_M < 0.2$')
p.ylabel(r"$\log_{10}$ (stellar mass / $M_\odot$ )")
#p.axvline(n.log10(3.), ls='dashed', label='z=0.1, 0.5, 1, 2')
#p.axvline(n.log10(2.), ls='dashed')#, label='z=1')
#p.axvline(n.log10(1.1),ls='dashed')#, label='z=0.1')
#p.axvline(n.log10(1.5),ls='dashed')#, label='z=0.5')
p.xlabel('redshift')
#p.xlabel(r'$\log_{10}(1+z)$')
#p.xscale('log')
p.legend(loc=0, frameon = False)
p.xlim((0.0, 1.4))
p.ylim((6.5, 12.5))
p.grid()
p.title('SDSS and eBOSS')
p.savefig(os.path.join(out_dir, "mass_redshift_mass_"+imf+"sdss_boss_02.jpg" ))
p.clf()


p.figure(1, (4.5, 4.5))
p.axes([0.2,0.2,0.7,0.7])
p.plot(zz_04, Ms_04, 'k,', rasterized=True, alpha=0.5, label=r' $\sigma_M < 0.4$')
#p.plot(zz_02, Ms_02, 'k,', rasterized=True, alpha=0.5, label=r' $\sigma_M < 0.2$')
p.ylabel(r"$\log_{10}$ (stellar mass / $M_\odot$ )")
#p.axvline(n.log10(3.), ls='dashed', label='z=0.1, 0.5, 1, 2')
#p.axvline(n.log10(2.), ls='dashed')#, label='z=1')
#p.axvline(n.log10(1.1),ls='dashed')#, label='z=0.1')
#p.axvline(n.log10(1.5),ls='dashed')#, label='z=0.5')
p.xlabel('redshift')
#p.xlabel(r'$\log_{10}(1+z)$')
#p.xscale('log')
p.legend(loc=0, frameon = False)
p.xlim((0.0, 1.4))
p.ylim((6.5, 12.5))
p.grid()
p.title('SDSS and eBOSS')
p.savefig(os.path.join(out_dir, "mass_redshift_mass_"+imf+"sdss_boss_04.jpg" ))
p.clf()

z_flg = 'ZQUALITY'
z_name = 'ZBEST'
deep2_zOk = (deep2[z_name] > 0.001) & (deep2[z_flg]>=2.) & (deep2[z_name] < 1.7) & (deep2['SSR']>0) & (deep2['TSR']>0) & (deep2['SSR']<=1.0001) & (deep2['TSR']<=1.0001)
deep2_sel_02 = (deep2_zOk) & (deep2[stellar_mass] < 10**14. ) & (deep2[stellar_mass] > 0. )  & (deep2[stellar_mass] >= deep2[stellar_mass+'_low'] ) & (deep2[stellar_mass] <= deep2[stellar_mass+'_up'] ) & ( - n.log10(deep2[stellar_mass+'_low'])  + n.log10(deep2[stellar_mass+'_up']) < 0.4 )
deep2_sel_04 = (deep2_zOk) & (deep2[stellar_mass] < 10**14. ) & (deep2[stellar_mass] > 0. )  & (deep2[stellar_mass] >= deep2[stellar_mass+'_low'] ) & (deep2[stellar_mass] <= deep2[stellar_mass+'_up'] ) & ( - n.log10(deep2[stellar_mass+'_low'])  + n.log10(deep2[stellar_mass+'_up']) < 0.8 )

Ms_02 = n.log10(deep2[stellar_mass][deep2_sel_02])
zz_02 = deep2[z_name][deep2_sel_02]
Ms_04 = n.log10(deep2[stellar_mass][deep2_sel_04])
zz_04 = deep2[z_name][deep2_sel_04]

w_deep2 = 1. / (deep2['TSR'] * deep2['SSR'])

print( "deep2 02",len(deep2_sel_02.nonzero()[0]))
print( "deep2 04",len(deep2_sel_04.nonzero()[0]))

p.figure(1, (4.5, 4.5))
p.axes([0.2,0.2,0.7,0.7])
#p.plot(zz_04, Ms_04, 'r+', rasterized=True, label=r' $\sigma_M < 0.4$', alpha=0.5)
p.plot(zz_02, Ms_02, 'k+', rasterized=True, label=r' $\sigma_M < 0.2$', alpha=0.5)
p.ylabel(r"$\log_{10}$ (stellar mass / $M_\odot$ )")
#p.axvline(n.log10(3.), ls='dashed', label='z=0.1, 0.5, 1, 2')
#p.axvline(n.log10(2.), ls='dashed')#, label='z=1')
#p.axvline(n.log10(1.1),ls='dashed')#, label='z=0.1')
#p.axvline(n.log10(1.5),ls='dashed')#, label='z=0.5')
p.xlabel('redshift')
#p.xlabel(r'$\log_{10}(1+z)$')
#p.xscale('log')
p.legend(loc=0, frameon = False)
p.xlim((0.0, 1.4))
p.ylim((6.5, 12.5))
p.grid()
p.title('DEEP2')
p.savefig(os.path.join(out_dir, "mass_redshift_mass_"+imf+"deep2_02.jpg" ))
p.clf()


p.figure(1, (4.5, 4.5))
p.axes([0.2,0.2,0.7,0.7])
p.plot(zz_04, Ms_04, 'k+', rasterized=True, label=r' $\sigma_M < 0.4$', alpha=0.5)
#p.plot(zz_02, Ms_02, 'kx', rasterized=True, label=r' $\sigma_M < 0.2$', alpha=0.5)
p.ylabel(r"$\log_{10}$ (stellar mass / $M_\odot$ )")
#p.axvline(n.log10(3.), ls='dashed', label='z=0.1, 0.5, 1, 2')
#p.axvline(n.log10(2.), ls='dashed')#, label='z=1')
#p.axvline(n.log10(1.1),ls='dashed')#, label='z=0.1')
#p.axvline(n.log10(1.5),ls='dashed')#, label='z=0.5')
p.xlabel('redshift')
#p.xlabel(r'$\log_{10}(1+z)$')
#p.xscale('log')
p.legend(loc=0, frameon = False)
p.xlim((0.0, 1.4))
p.ylim((6.5, 12.5))
p.grid()
p.title('DEEP2')
p.savefig(os.path.join(out_dir, "mass_redshift_mass_"+imf+"deep2_04.jpg" ))
p.clf()

z_bins = n.arange(0.,1.4,0.05)
p.figure(2, (6.5, 3.5))
p.axes([0.12,0.18,0.8,0.73])
p.hist(zz_02_boss, bins = z_bins, histtype='step', label=r'eBOSS', ls='solid', color='r')
p.hist(zz_04_boss, bins = z_bins, histtype='step', ls='dashed', color='r')
p.hist(zz_02_sdss, bins = z_bins, histtype='step', label='SDSS', ls='solid', color='k')
p.hist(zz_04_sdss, bins = z_bins, histtype='step', ls='dashed', color='k')
p.hist(zz_02     , bins = z_bins, histtype='step', label='DEEP2', ls='solid', color='b')
p.hist(zz_04     , bins = z_bins, histtype='step', ls='dashed', color='b')
p.ylabel(r"N(dz=0.05)")
p.xlabel('redshift')
p.yscale('log')
p.legend(loc=0, frameon = False)
p.xlim((0.0, 1.4))
#p.ylim((1, 12.5))
p.grid()
p.savefig(os.path.join(out_dir, "redshift_distribution.jpg" ))
p.clf()

import sys
sys.exit()

def plot_all_prognames(hdus=hdus, imf=imf, prefix=prefix, out_dir = out_dir, redshift_reliable=redshift_reliable ) :
	stellar_mass = imf+'_stellar_mass'
	out_dir = os.path.join(os.environ['DATA_DIR'], 'spm', 'results', 'catalogs', imf)

	mass_reliable = (hdus[1].data[stellar_mass] > 0 ) & ( hdus[1].data[stellar_mass] < 13. ) & ( abs(hdus[1].data[stellar_mass + '_err_plus'] - hdus[1].data[stellar_mass]) < 0.4 ) & ( abs(hdus[1].data[stellar_mass + '_err_minus'] - hdus[1].data[stellar_mass]) < 0.4 )

	#good_plates = (hdus[1].data['PLATEQUALITY']=='good') &(hdus[1].data['TARGETTYPE']=='science')

	all_names = set(hdus[1].data['PROGRAMNAME'])
	all_names_arr = n.array(list(all_names))

	for ii in range(len(all_names_arr)):
		selection = (mass_reliable) & (redshift_reliable) & (hdus[1].data['PROGRAMNAME']==all_names_arr[ii])
		N_occ = len(selection.nonzero()[0])
		print all_names_arr[ii], N_occ
		if N_occ>1:
			p.figure(1, (4.5, 4.5))
			p.axes([0.2,0.2,0.7,0.7])
			p.plot(n.log10(1.+hdus[1].data['Z'][selection]), hdus[1].data[stellar_mass][selection], 'k+', rasterized=True, alpha=0.5) #, label=all_names_arr[ii]
			p.ylabel(r'$\log_{10}$ (stellar mass '+imf+r" / $M_\odot$ )")
			p.axvline(n.log10(3.), ls='dashed', label='z=0.1, 0.5, 1, 2')
			p.axvline(n.log10(2.), ls='dashed')#, label='z=1')
			p.axvline(n.log10(1.1),ls='dashed')#, label='z=0.1')
			p.axvline(n.log10(1.5),ls='dashed')#, label='z=0.5')
			p.xlabel(r'$\log_{10}(1+z)$')
			#p.xscale('log')
			p.legend(loc=0, frameon = False)
			p.xlim((0.0, 0.7))
			p.ylim((6.5, 12.5))
			p.grid()
			p.title('N='+str(N_occ))
			p.savefig(os.path.join(out_dir, prefix+"_"+all_names_arr[ii]+"_redshift_mass_"+imf+".jpg" ))
			p.clf()

			p.figure(2, (4.5, 4.5))
			p.axes([0.2,0.2,0.7,0.7])
			#p.subplot(111, projection="mollweide")
			#p.plot((hdus[1].data['PLUG_RA'][selection]-180.)*n.pi/180., hdus[1].data['PLUG_DEC'][selection]*n.pi/180., 'k+', rasterized=True) # , label=all_names_arr[ii]
			p.plot(hdus[1].data['PLUG_RA'][selection], hdus[1].data['PLUG_DEC'][selection], 'k+', rasterized=True) # , label=all_names_arr[ii]
			p.title(all_names_arr[ii])#+', Ngal='+str(N_occ))
			p.xlim((0.0, 360.))
			p.ylim((-20., 85.))
			p.xlabel(r'ra [deg]')
			p.ylabel(r'dec [deg]')
			#p.legend(loc=0, frameon = False)
			p.grid()
			p.savefig(os.path.join(out_dir, prefix+"_"+all_names_arr[ii]+"_ra_dec_"+imf+".jpg" ))
			p.clf()



plot_all_prognames(hdus=hdus, imf=imf, prefix=prefix, redshift_reliable=redshift_reliable )
"""
imf = 'Salpeter'
plot_all_prognames(hdus=hdus, imf=imf, prefix=prefix,  redshift_reliable=redshift_reliable )

imf = 'Kroupa'
plot_all_prognames(hdus=hdus, imf=imf, prefix=prefix, redshift_reliable=redshift_reliable )
"""
prefix = 'SDSS'
hdus = fits.open(os.path.join(os.environ['DATA_DIR'], 'spm', 'firefly', 'FireflyGalaxySdss26.fits'))
redshift_reliable = (hdus[1].data['SN_MEDIAN_ALL'] > 0.1 ) & (hdus[1].data['CLASS'] == "GALAXY") & (hdus[1].data['Z'] >= 0) & ( hdus[1].data['Z_ERR'] >= 0) & (hdus[1].data['ZWARNING'] == 0) & (hdus[1].data['Z'] > hdus[1].data['Z_ERR'] )

imf = 'Chabrier'
plot_all_prognames(hdus=hdus, imf=imf, prefix=prefix, redshift_reliable=redshift_reliable )
"""
imf = 'Salpeter'
plot_all_prognames(hdus=hdus, imf=imf, prefix=prefix, redshift_reliable=redshift_reliable )

imf = 'Kroupa'
plot_all_prognames(hdus=hdus, imf=imf, prefix=prefix, redshift_reliable=redshift_reliable )
"""




import astropy.cosmology as co
cosmo=co.Planck15
import astropy.io.fits as fits
import matplotlib.pyplot as p
import numpy as n
import os
import sys
from scipy.integrate import quad
# stat functions
ld = lambda selection : len(selection.nonzero()[0])

# global cosmo quantities
z_min = 0.2
z_max = 0.5
volume_per_deg2 = ( cosmo.comoving_volume(z_max) -  cosmo.comoving_volume(z_min) ) * n.pi / 129600.
volume_per_deg2_val = volume_per_deg2.value

smf_ilbert13 = lambda M, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s : ( phi_1s * (M/M_star) ** alpha_1s + phi_2s * (M/M_star) ** alpha_2s ) * n.e ** (-M/M_star) * (M/ M_star)
path_ilbert13_SMF = os.path.join(os.environ['OBS_REPO'], 'spm/literature', "ilbert_2013_mass_function_params.txt")
zmin, zmax, N, M_comp, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s, log_rho_s = n.loadtxt(path_ilbert13_SMF, unpack=True)
smf01 = lambda mass : smf_ilbert13( mass , 10**M_star[0], phi_1s[0]*10**(-3), alpha_1s[0], phi_2s[0]*10**(-3), alpha_2s[0] )
#print 10**M_star[0], phi_1s[0]*10**(-3), alpha_1s[0], phi_2s[0]*10**(-3), alpha_2s[0]

fun = lambda mass :smf01(mass) #  mass * 
masses = 10**n.arange(7,12,0.1)
total_mass_per_unit_volume = n.array([quad(fun, mmin, 10**13)[0] for mmin in masses])

volume_sdss = 7900.   * volume_per_deg2_val 
volume_boss = 10000.  * volume_per_deg2_val
volume_cosmos = 1.52  * volume_per_deg2_val
volume_deep2 = 2.0  * volume_per_deg2_val


out_dir = os.path.join(os.environ['OBS_REPO'], 'spm', 'results', 'mass-density')

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

# OPENS THE CATALOGS
deep2   = fits.open(path_2_deep2_cat)[1].data
sdss   = fits.open(path_2_sdss_cat)[1].data
boss   = fits.open(path_2_eboss_cat)[1].data
cosmos = fits.open(path_2_cosmos_cat)[1].data


imf = imfs[0]
stellar_mass = imf+'stellar_mass'

redshift_reliable_boss =  (boss['CLASS_NOQSO'] == "GALAXY") & ( boss['Z_ERR_NOQSO'] > 0.0) & (boss['ZWARNING_NOQSO'] == 0) & (boss['Z_NOQSO']>0.001) & (boss['Z_NOQSO'] > boss['Z_ERR_NOQSO'] ) # (boss['SN_MEDIAN_ALL'] > 0.1 ) & 
error_reliable_boss = (boss[stellar_mass+'_up'] > boss[stellar_mass+'_low'] ) & (boss[stellar_mass+'_up'] > 0. ) & ( boss[stellar_mass+'_low'] > 0. ) & (boss[stellar_mass+'_up'] < 1e14 ) & ( boss[stellar_mass+'_low'] < 1e14 ) 
mass_reliable_boss_02 = (boss[stellar_mass] > 1e6 ) & ( boss[stellar_mass] < 1e14 ) & ((n.log10(boss[stellar_mass+'_up']) - n.log10(boss[stellar_mass+'_low']))/2. < 0.2 )
mass_reliable_boss_04 = (boss[stellar_mass] > 1e6 ) & ( boss[stellar_mass] < 1e14 ) & ((n.log10(boss[stellar_mass+'_up']) - n.log10(boss[stellar_mass+'_low']))/2. < 0.4 )
ok_boss_02 = (error_reliable_boss) & (mass_reliable_boss_02) & (redshift_reliable_boss)& (boss['Z_NOQSO']>z_min)& (boss['Z_NOQSO']<z_max)
ok_boss_04 = (error_reliable_boss) & (mass_reliable_boss_04) & (redshift_reliable_boss)& (boss['Z_NOQSO']>z_min)& (boss['Z_NOQSO']<z_max)

redshift_reliable_sdss =  (sdss['CLASS'] == "GALAXY")       & ( sdss['Z_ERR'] > 0.0)       & (sdss['ZWARNING'] == 0)       & (sdss['Z'] > 0.001) & (sdss['Z'] > sdss['Z_ERR'] ) # (sdss['SN_MEDIAN_ALL'] > 0.1 ) &
error_reliable_sdss = (sdss[stellar_mass+'_up'] > sdss[stellar_mass+'_low'] ) & (sdss[stellar_mass+'_up'] > 0. ) & ( sdss[stellar_mass+'_low'] > 0. ) & (sdss[stellar_mass+'_up'] < 1e14 ) & ( sdss[stellar_mass+'_low'] < 1e14 ) 
mass_reliable_sdss_02 = (sdss[stellar_mass] > 1e6 ) & ( sdss[stellar_mass] < 1e14 ) & ((n.log10(sdss[stellar_mass+'_up']) - n.log10(sdss[stellar_mass+'_low']))/2. < 0.2 )
mass_reliable_sdss_04 = (sdss[stellar_mass] > 1e6 ) & ( sdss[stellar_mass] < 1e14 ) & ((n.log10(sdss[stellar_mass+'_up']) - n.log10(sdss[stellar_mass+'_low']))/2. < 0.4 )
ok_sdss_02 = (error_reliable_sdss) & (mass_reliable_sdss_02) & (redshift_reliable_sdss)& (sdss['Z']>z_min)& (sdss['Z']<z_max)
ok_sdss_04 = (error_reliable_sdss) & (mass_reliable_sdss_04) & (redshift_reliable_sdss)& (sdss['Z']>z_min)& (sdss['Z']<z_max)

z_flg = 'ZQUALITY'
z_name = 'ZBEST'
deep2_zOk = (deep2[z_name] > 0.001) & (deep2[z_flg]>=2.) & (deep2[z_name] < 1.7) & (deep2['SSR']>0) & (deep2['TSR']>0) & (deep2['SSR']<=1.0001) & (deep2['TSR']<=1.0001)
ok_deep2_02 = (deep2_zOk) & (deep2[stellar_mass] < 10**14. ) & (deep2[stellar_mass] > 0. )  & (deep2[stellar_mass] >= deep2[stellar_mass+'_low'] ) & (deep2[stellar_mass] <= deep2[stellar_mass+'_up'] ) & ( - n.log10(deep2[stellar_mass+'_low'])  + n.log10(deep2[stellar_mass+'_up']) < 0.4 )
ok_deep2_04 = (deep2_zOk) & (deep2[stellar_mass] < 10**14. ) & (deep2[stellar_mass] > 0. )  & (deep2[stellar_mass] >= deep2[stellar_mass+'_low'] ) & (deep2[stellar_mass] <= deep2[stellar_mass+'_up'] ) & ( - n.log10(deep2[stellar_mass+'_low'])  + n.log10(deep2[stellar_mass+'_up']) < 0.8 )

total_mass_boss_02 = n.array([ n.sum(boss[stellar_mass][ok_boss_02 & (boss[stellar_mass] > mmin)]) for mmin in masses ])
total_mass_boss_04 = n.array([ n.sum(boss[stellar_mass][ok_boss_04 & (boss[stellar_mass] > mmin)]) for mmin in masses ])

total_mass_sdss_02 = n.array([ n.sum(sdss[stellar_mass][ok_sdss_02 & (sdss[stellar_mass] > mmin)]) for mmin in masses ])
total_mass_sdss_04 = n.array([ n.sum(sdss[stellar_mass][ok_sdss_04 & (sdss[stellar_mass] > mmin)]) for mmin in masses ])

total_mass_deep2_02 = n.array([ n.sum(deep2[stellar_mass][ok_deep2_02 & (deep2[stellar_mass] > mmin)]) for mmin in masses ])
total_mass_deep2_04 = n.array([ n.sum(deep2[stellar_mass][ok_deep2_04 & (deep2[stellar_mass] > mmin)]) for mmin in masses ])


p.figure(1, (4.5, 4.5))
p.axes([0.2,0.2,0.7,0.7])
total = total_mass_per_unit_volume * volume_boss
p.plot(masses, total_mass_boss_02/total, label='BOSS 0.2' )
#p.plot(masses, total_mass_boss_04/total, label='BOSS 0.4' )
total = total_mass_per_unit_volume * volume_sdss
p.plot(masses, total_mass_sdss_02/total, label='SDSS 0.2' )
#p.plot(masses, total_mass_sdss_04/total, label='SDSS 0.4' )
total = total_mass_per_unit_volume * volume_deep2
p.plot(masses, total_mass_deep2_02/total, label='DEEP2 0.2' )
#p.plot(masses, total_mass_deep2_04/total, label='DEEP2 0.4' )
p.legend(frameon=False, loc=0)
p.xlabel('Mmin')
p.ylabel('Mass fraction (M>Mmin)')
p.title(str(z_min)+'<z<'+str(z_max))
p.xscale('log')
p.yscale('log')
p.ylim((0.001,2))
p.xlim((1e7,10**(12.5)))
p.grid()
p.savefig(os.path.join(out_dir, "mass_density"+imf+"_"+str(z_min)+'_z_'+str(z_max)+".jpg" ))
p.show()
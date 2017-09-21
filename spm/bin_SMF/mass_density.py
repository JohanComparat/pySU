import astropy.cosmology as co
cosmo=co.Planck15
import astropy.io.fits as fits
import matplotlib
matplotlib.rcParams['agg.path.chunksize'] = 2000000
matplotlib.rcParams.update({'font.size': 12})
import matplotlib.pyplot as p
import numpy as n
import os
import sys
from scipy.integrate import quad
# stat functions
ld = lambda selection : len(selection.nonzero()[0])

# global cosmo quantities
#run mass_density.py 0.2 0.5
#run mass_density.py 0.5 0.8
#run mass_density.py 0.8 1.1
#run mass_density.py 1.1 1.5
z_min = float(sys.argv[1])
z_max = float(sys.argv[2])

volume_per_deg2 = ( cosmo.comoving_volume(z_max) -  cosmo.comoving_volume(z_min) ) * n.pi / 129600.
volume_per_deg2_val = volume_per_deg2.value

smf_ilbert13 = lambda M, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s : ( phi_1s * (M/M_star) ** alpha_1s + phi_2s * (M/M_star) ** alpha_2s ) * n.e ** (-M/M_star) * (M/ M_star)
path_ilbert13_SMF = os.path.join(os.environ['OBS_REPO'], 'spm/literature', "ilbert_2013_mass_function_params.txt")
zmin, zmax, N, M_comp, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s, log_rho_s = n.loadtxt(path_ilbert13_SMF, unpack=True)

if z_min == 0.2 :
  smf01 = lambda mass : smf_ilbert13( mass , 10**M_star[0], phi_1s[0]*10**(-3), alpha_1s[0], phi_2s[0]*10**(-3), alpha_2s[0] )
  area_deep2 = 0.5
if z_min == 0.5 :
  smf01 = lambda mass : smf_ilbert13( mass , 10**M_star[1], phi_1s[1]*10**(-3), alpha_1s[1], phi_2s[1]*10**(-3), alpha_2s[1] )
  area_deep2 = 2.0 
if z_min == 0.8 :
  smf01 = lambda mass : smf_ilbert13( mass , 10**M_star[2], phi_1s[2]*10**(-3), alpha_1s[2], phi_2s[2]*10**(-3), alpha_2s[2] )
  area_deep2 = 2.0 
if z_min == 1.1 :
  smf01 = lambda mass : smf_ilbert13( mass , 10**M_star[3], phi_1s[3]*10**(-3), alpha_1s[3], phi_2s[3]*10**(-3), alpha_2s[3] )
  area_deep2 = 2.0 
  
#print 10**M_star[0], phi_1s[0]*10**(-3), alpha_1s[0], phi_2s[0]*10**(-3), alpha_2s[0]

fun = lambda mass :smf01(mass) #  mass * 
masses = 10**n.arange(7,12,0.1)
total_mass_per_unit_volume = n.array([quad(fun, mmin, 10**13)[0] for mmin in masses])


area_sdss = 7900. 
area_boss = 10000.
area_cosmos = 1.52
 

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


def get_basic_stat_DR14(catalog, z_name, z_err_name, class_name, zwarning,  zflg_val, prefix, err_max):
    catalog_zOk =(catalog[z_err_name] > 0.) & (catalog[z_name] > catalog[z_err_name])  & (catalog[class_name]=='GALAXY')  & (catalog[zwarning]==zflg_val)
    catalog_stat = (catalog_zOk) & (catalog[z_name] > z_min) & (catalog[z_name] < z_max) 
    catalog_sel = (catalog_stat) & (catalog[prefix+'stellar_mass'] < 10**14. ) & (catalog[prefix+'stellar_mass'] > 0 )  & (catalog[prefix+'stellar_mass'] > catalog[prefix+'stellar_mass_low'] ) & (catalog[prefix+'stellar_mass'] < catalog[prefix+'stellar_mass_up'] ) & ( - n.log10(catalog[prefix+'stellar_mass_low'])  + n.log10(catalog[prefix+'stellar_mass_up']) < err_max )
    m_catalog = n.log10(catalog[prefix+'stellar_mass'])
    w_catalog =  n.ones_like(catalog[prefix+'stellar_mass'])
    return catalog_sel, m_catalog, w_catalog

def get_basic_stat_DEEP2(deep2, IMF, err_max):
    z_flg = 'ZQUALITY'
    z_name = 'ZBEST'
    stellar_mass = IMF+'stellar_mass'
    deep2_zOk = (deep2[z_name] > z_min) & (deep2[z_name] < z_max) & (deep2[z_flg]>=2.) & (deep2['SSR']>0) & (deep2['TSR']>0) & (deep2['SSR']<=1.0001) & (deep2['TSR']<=1.0001)
    ok_deep2_02 = (deep2_zOk) & (deep2[stellar_mass] < 10**14. ) & (deep2[stellar_mass] > 0. ) & ( - n.log10(deep2[stellar_mass+'_low'])  + n.log10(deep2[stellar_mass+'_up']) < err_max)
    return ok_deep2_02, n.log10(deep2[stellar_mass][ok_deep2_02]), 1./(deep2['SSR'][ok_deep2_02]*deep2['TSR'][ok_deep2_02])

    
def get_hist(masses, weights, mbins):
    NN = n.histogram(masses, mbins)[0]
    NW = n.histogram(masses, mbins, weights = weights)[0]
    xx = (mbins[1:] + mbins[:-1])/2.
    return xx, NW, NN**(-0.5)*NW


dlog10m = 0.25
mbins = n.arange(8,12.5,dlog10m)

def plot_smf_b(IMF="Chabrier_ELODIE_", err_max=0.4):
	boss_sel, boss_m, boss_w = get_basic_stat_DR14(boss, 'Z_NOQSO', 'Z_ERR_NOQSO', 'CLASS_NOQSO', 'ZWARNING_NOQSO', 0., IMF, err_max)
	x, y, ye = get_hist(boss_m[boss_sel], weights = boss_w[boss_sel]/(dlog10m*n.log(10)*area_boss*volume_per_deg2_val), mbins = mbins)
	return x, y, ye
	#p.errorbar(x, y, yerr = ye, label=IMF[:-1], lw=1)

def plot_smf_s(IMF="Chabrier_ELODIE_", err_max=0.4):
	boss_sel, boss_m, boss_w = get_basic_stat_DR14(sdss, 'Z', 'Z_ERR', 'CLASS', 'ZWARNING', 0., IMF, err_max)
	x, y, ye = get_hist(boss_m[boss_sel], weights = boss_w[boss_sel]/(dlog10m*n.log(10)*area_sdss*volume_per_deg2_val), mbins = mbins)
	return x, y, ye
	#p.errorbar(x, y, yerr = ye, label=IMF[:-1], lw=1)

def plot_smf_d(IMF="Chabrier_ELODIE_", err_max=0.4, area_deep2=0.5):
	boss_sel, boss_m, boss_w = get_basic_stat_DEEP2(deep2, IMF, err_max)
	x, y, ye = get_hist(boss_m, weights = boss_w/(dlog10m*n.log(10)*area_deep2*volume_per_deg2_val), mbins = mbins)
	return x, y, ye


xa, ya, yea = plot_smf_b("Chabrier_ELODIE_", 0.2*2)
xb, yb, yeb = plot_smf_b("Chabrier_MILES_", 0.2*2)
xc, yc, yec = plot_smf_b("Chabrier_STELIB_", 0.2*2)

xd, yd, yed = plot_smf_s("Chabrier_ELODIE_", 0.2*2)
xe, ye, yee = plot_smf_s("Chabrier_MILES_", 0.2*2)
xf, yf, yef = plot_smf_s("Chabrier_STELIB_", 0.2*2)

xg, yg, yeg = plot_smf_d("Chabrier_ELODIE_", 0.2*2, area_deep2=area_deep2)
xh, yh, yeh = plot_smf_d("Chabrier_MILES_", 0.2*2, area_deep2=area_deep2)
xi, yi, yei = plot_smf_d("Chabrier_STELIB_", 0.2*2, area_deep2=area_deep2)

p.figure(1, (4.5,4.5))
p.axes([0.19,0.17,0.74,0.72])
p.fill_between( mbins, y1=smf01(10**mbins)*0.77, y2=smf01(10**mbins)*1.23, color='g', alpha=0.5)
p.plot(mbins, smf01(10**mbins), label='Ilbert 13', color='g')

p.fill_between( xa, y1=n.min([ya-yea, yb-yeb, yc-yec], axis=0), y2=n.max([ya+yea, yb+yeb, yc+yec], axis=0), color='r', alpha=0.5)
p.plot(xa, n.mean([ya, yb, yc], axis=0), label=r'BOSS, eBOSS', color='r')

p.fill_between( xa, y1=n.min([yd-yed, ye-yee, yf-yef], axis=0), y2=n.max([yd+yed, ye+yee, yf+yef], axis=0), color='b', alpha=0.5)
p.plot(xa, n.mean([yd, ye, yf], axis=0), label=r'SDSS', color='b')

p.fill_between( xa, y1=n.min([yg-yeg, yh-yeh, yi-yei], axis=0), y2=n.max([yg+yeg, yh+yeh, yi+yei], axis=0), color='k', alpha=0.5)
p.plot(xa, n.mean([yg, yh, yi], axis=0), label=r'DEEP2', color='k')
p.title('Chabrier IMF '+str(z_min)+'<z<'+str(z_max))
p.xlabel(r"$\log_{10}$ (M / $M_\odot$ )")
p.ylabel(r'$\Phi(M)$ [Mpc$^{-3}$ dex$^{-1}$]')
p.yscale('log')
p.legend(loc=6, frameon = False)
p.ylim((1e-9, 1e-2))
p.xlim((8.5, 12.5))
p.grid()
p.savefig(os.path.join(out_dir, "firefly_SMF_BOSS_"+str(z_min)+'_z_'+str(z_max)+".jpg" ))
p.clf()

sys.exit()


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
p.clf()

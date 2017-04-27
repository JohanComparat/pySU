import astropy.cosmology as co
aa=co.Planck15
import astropy.io.fits as fits
import matplotlib.pyplot as p
import numpy as n
import os

import sys

# global cosmo quantities
z_min = 0. # float(sys.argv[1])
z_max = 6. # float(sys.argv[2])
#imf = 'kroupa'
#lO2_min = float(sys.argv[3]) # 'salpeter'

SNlimit = 5

smf_ilbert13 = lambda M, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s : ( phi_1s * (M/M_star) ** alpha_1s + phi_2s * (M/M_star) ** alpha_2s ) * n.e ** (-M/M_star) * (M/ M_star)

ff_dir = os.path.join(os.environ['DATA_DIR'], 'spm', 'firefly')
ll_dir = os.path.join(os.environ['DATA_DIR'], 'spm', 'literature')
co_dir = os.path.join(os.environ['DATA_DIR'], 'COSMOS' )
sdss_dir = os.path.join(os.environ['DATA_DIR'], 'SDSS')

out_dir = os.path.join(os.environ['DATA_DIR'], 'spm', 'results')

path_2_cosmos_cat = os.path.join( co_dir, "photoz_vers2.0_010312.fits")

#path_2_sdss_cat = os.path.join( ff_dir, "FireflyGalaxySdssDR14.fits" )
#path_2_eboss_cat = os.path.join( ff_dir, "FireflyGalaxyEbossDR14.fits" )

path_2_spall_sdss_dr12_cat = os.path.join( sdss_dir, "specObj-SDSS-dr12.fits" )
path_2_spall_sdss_dr14_cat = os.path.join( sdss_dir, "specObj-SDSS-dr14.fits" )
path_2_spall_boss_dr12_cat = os.path.join( sdss_dir, "specObj-BOSS-dr12.fits" )
path_2_spall_boss_dr14_cat = os.path.join( sdss_dir, "specObj-BOSS-dr14.fits" )

#print "SDSS spAll DR12", len(fits.open(path_2_spall_sdss_dr12_cat)[1].data)
print "SDSS spAll DR14", len(fits.open(path_2_spall_sdss_dr14_cat)[1].data)
print "BOSS spAll DR12",len(fits.open(path_2_spall_boss_dr12_cat)[1].data)
print "BOSS spAll DR14",len(fits.open(path_2_spall_boss_dr14_cat)[1].data)

path_2_sdss_kroupa_cat = os.path.join( ff_dir, "spAll-SDSS-FireflyKroupa.fits" )
path_2_sdss_salpeter_cat = os.path.join( ff_dir, "spAll-SDSS-FireflySalpeter.fits" )
path_2_eboss_kroupa_cat = os.path.join( ff_dir, "spAll-eBOSS-FireflyKroupa.fits" )
path_2_eboss_salpeter_cat = os.path.join( ff_dir, "spAll-eBOSS-FireflySalpeter.fits" )

path_2_pS_salpeter_cat = os.path.join( ll_dir, "portsmouth_stellarmass_starforming_salp-26.fits.gz" )
path_2_pB_salpeter_cat = os.path.join( ll_dir, "portsmouth_stellarmass_starforming_salp-DR12-boss.fits.gz" )

path_2_pS_kroupa_cat = os.path.join( ll_dir, "portsmouth_stellarmass_starforming_krou-26.fits.gz" )
path_2_pB_kroupa_cat = os.path.join( ll_dir, "portsmouth_stellarmass_starforming_krou-DR12-boss.fits.gz" )

path_2_ppS_kroupa_cat = os.path.join( ll_dir, "portsmouth_stellarmass_passive_krou-26.fits")
path_2_ppB_kroupa_cat = os.path.join( ll_dir, "portsmouth_stellarmass_passive_krou-DR12.fits")

path_2_F16_cat = os.path.join( sdss_dir, "RA_DEC_z_w_fluxOII_Mstar_grcol_Mr_lumOII.dat" )

RA, DEC, z, weigth, O2flux, M_star, gr_color, Mr_5logh, O2luminosity = n.loadtxt(path_2_F16_cat, unpack=True)

cosmos = fits.open(path_2_cosmos_cat)[1].data
sdss_k   = fits.open(path_2_sdss_kroupa_cat)[1].data
boss_k   = fits.open(path_2_eboss_kroupa_cat)[1].data
boss_s   = fits.open(path_2_eboss_salpeter_cat)[1].data
sdss_s   = fits.open(path_2_sdss_salpeter_cat)[1].data

sdss_12_portSF_kr   = fits.open(path_2_pS_kroupa_cat)[1].data
boss_12_portSF_kr   = fits.open(path_2_pB_kroupa_cat)[1].data

sdss_12_portPA_kr   = fits.open(path_2_ppS_kroupa_cat)[1].data
boss_12_portPA_kr   = fits.open(path_2_ppB_kroupa_cat)[1].data

sdss_12_portSF_sa   = fits.open(path_2_pS_salpeter_cat)[1].data
boss_12_portSF_sa   = fits.open(path_2_pB_salpeter_cat)[1].data


path_ilbert13_SMF = os.path.join(ll_dir, "ilbert_2013_mass_function_params.txt")
zmin, zmax, N, M_comp, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s, log_rho_s = n.loadtxt(os.path.join( ll_dir, "ilbert_2013_mass_function_params.txt"), unpack=True)

#smfs_ilbert13 = n.array([lambda mass : smf_ilbert13( mass , 10**M_star[ii], phi_1s[ii]*10**(-3), alpha_1s[ii], phi_2s[ii]*10**(-3), alpha_2s[ii] ) for ii in range(len(M_star)) ])

smf01 = lambda mass : smf_ilbert13( mass , 10**M_star[0], phi_1s[0]*10**(-3), alpha_1s[0], phi_2s[0]*10**(-3), alpha_2s[0] )
print 10**M_star[0], phi_1s[0]*10**(-3), alpha_1s[0], phi_2s[0]*10**(-3), alpha_2s[0]

volume_per_deg2 = ( aa.comoving_volume(z_max) -  aa.comoving_volume(z_min) ) * n.pi / 129600.
volume_per_deg2_val = volume_per_deg2.value

# global spm quantities

# stat functions
ld = lambda selection : len(selection.nonzero()[0])

area_sdss14 = 7900.    
area_sdss12 = 7900.

area_boss14 = 10000.
area_boss12 = 10000.

area_cosmos = 1.52


def get_basic_stat_anyCat(catalog_name, z_name, z_err_name, name, zflg_val):
    catalog = fits.open(catalog_name)[1].data
    catalog_zOk =(catalog[z_err_name] > 0.) & (catalog[z_name] > catalog[z_err_name]) 
    catalog_stat = (catalog_zOk) #& (catalog[z_name] > z_min) & (catalog[z_name] < z_max) 
    catalog_sel = (catalog_stat) & (catalog['LOGMASS'] < 14. ) & (catalog['LOGMASS'] > 0 ) & (catalog['MAXLOGMASS'] - catalog['MINLOGMASS'] <0.4) & (catalog['LOGMASS'] < catalog['MAXLOGMASS'] ) & (catalog['LOGMASS'] > catalog['MINLOGMASS'] )
    m_catalog = catalog['LOGMASS']
    w_catalog =  n.ones_like(catalog['LOGMASS'])
    print catalog_name, "& - & $", ld(catalog_zOk),"$ & $", ld(catalog_sel),"$ \\\\"
    #return catalog_sel, m_catalog, w_catalog

def get_basic_stat_DR12(catalog, z_name, z_err_name, name, zflg_val):
    catalog_zOk =(catalog[z_err_name] > 0.) & (catalog[z_name] > catalog[z_err_name]) 
    catalog_stat = (catalog_zOk) #& (catalog[z_name] > z_min) & (catalog[z_name] < z_max) 
    catalog_sel = (catalog_stat) & (catalog['LOGMASS'] < 14. ) & (catalog['LOGMASS'] > 0 ) & (catalog['MAXLOGMASS'] - catalog['MINLOGMASS'] <0.4) & (catalog['LOGMASS'] < catalog['MAXLOGMASS'] ) & (catalog['LOGMASS'] > catalog['MINLOGMASS'] )
    m_catalog = catalog['LOGMASS']
    w_catalog =  n.ones_like(catalog['LOGMASS'])
    print name, "& - & $", ld(catalog_zOk),"$ & $", ld(catalog_sel),"$ \\\\"
    return catalog_sel, m_catalog, w_catalog

def get_basic_stat_DR14(catalog, z_name, z_err_name, class_name, zwarning, name, zflg_val):
    catalog_zOk =(catalog[z_err_name] > 0.) & (catalog[z_name] > catalog[z_err_name])  & (catalog[class_name]=='GALAXY')  & (catalog[zwarning]==zflg_val)
    catalog_stat = (catalog_zOk) #& (catalog[z_name] > z_min) & (catalog[z_name] < z_max) 
    catalog_sel = (catalog_stat) & (catalog['stellar_mass'] < 14. ) & (catalog['stellar_mass'] > 0 )  & (catalog['stellar_mass'] > catalog['stellar_mass_err_plus'] ) & (catalog['stellar_mass'] > catalog['stellar_mass_err_minus'] ) & (catalog['stellar_mass_err_minus']  + catalog['stellar_mass_err_plus'] <0.4)
    m_catalog = catalog['stellar_mass']
    w_catalog =  n.ones_like(catalog['stellar_mass'])
    print name, '& $',len(catalog), "$ & $", ld(catalog_zOk),"$ & $", ld(catalog_sel),"$ \\\\"
    return catalog_sel, m_catalog, w_catalog

def get_basic_stat_FF(catalog, z_name, z_err_name, z_flg, name, zflg_val, imf='kroupa'):
    catalog_zOk =(catalog[z_err_name] > 0.) & (catalog[z_name] > catalog[z_err_name]) & (catalog[z_flg]==zflg_val) 
    catalog_stat = (catalog_zOk) & (catalog[z_name] > z_min) & (catalog[z_name] < z_max) 
    catalog_sel = (catalog_stat) & (catalog['stellar_mass_'+imf] < 14. ) & (catalog['stellar_mass_'+imf] > 0 )  & (catalog['stellar_mass_'+imf] > catalog['stellar_mass_err_plus_'+imf] ) & (catalog['stellar_mass_'+imf] > catalog['stellar_mass_err_minus_'+imf] ) & (catalog['stellar_mass_err_minus_'+imf]  + catalog['stellar_mass_err_plus_'+imf] <0.4)
    m_catalog = catalog['stellar_mass_'+imf]
    w_catalog =  n.ones_like(catalog['stellar_mass_'+imf])
    print name, '& $',len(catalog), "$ & $", ld(catalog_zOk),"$ & $", ld(catalog_stat), "\\;(", ld(catalog_sel),")$ \\\\"
    return catalog_sel, m_catalog, w_catalog
    
def get_hist(masses, weights, mbins):
    NN = n.histogram(masses, mbins)[0]
    NW = n.histogram(masses, mbins, weights = weights)[0]
    xx = (mbins[1:] + mbins[:-1])/2.
    return xx, NW, NN**(-0.5)*NW


boss12_sel, boss12_m, boss12_w = get_basic_stat_DR14(boss_k, 'Z_NOQSO', 'Z_ERR_NOQSO', 'CLASS_NOQSO', 'ZWARNING_NOQSO', 'FIREFLY Kroupa & BOSS & 14 ', 0.)
boss12_sel, boss12_m, boss12_w = get_basic_stat_DR14(boss_s, 'Z_NOQSO', 'Z_ERR_NOQSO', 'CLASS_NOQSO', 'ZWARNING_NOQSO', 'FIREFLY Salpeter & BOSS & 14', 0.)
sdss12_sel, sdss12_m, sdss12_w = get_basic_stat_DR12(boss_12_portSF_kr, 'Z', 'Z_ERR', 'Portsmouth SF Kroupa   & BOSS & 12 ', 0.)
sdss12_sel, sdss12_m, sdss12_w = get_basic_stat_DR12(boss_12_portPA_kr, 'Z', 'Z_ERR', 'Portsmouth Passive Kroupa   & BOSS & 12 ', 0.)
sdss12_sel, sdss12_m, sdss12_w = get_basic_stat_DR12(boss_12_portSF_sa, 'Z', 'Z_ERR', 'Portsmouth SF Salpeter & BOSS & 12 ', 0.)

sdss12_sel, sdss12_m, sdss12_w = get_basic_stat_DR14(sdss_k, 'Z', 'Z_ERR', 'CLASS', 'ZWARNING', 'FIREFLY Kroupa & SDSS & 14 ', 0.)
sdss12_sel, sdss12_m, sdss12_w = get_basic_stat_DR14(sdss_s, 'Z', 'Z_ERR', 'CLASS', 'ZWARNING', 'FIREFLY Salpeter & SDSS & 14', 0.)
sdss12_sel, sdss12_m, sdss12_w = get_basic_stat_DR12(sdss_12_portSF_kr, 'Z', 'Z_ERR', 'Portsmouth SF Kroupa   & SDSS & 12 ', 0.)
sdss12_sel, sdss12_m, sdss12_w = get_basic_stat_DR12(sdss_12_portPA_kr, 'Z', 'Z_ERR', 'Portsmouth Passive Kroupa   & SDSS & 12 ', 0.)
sdss12_sel, sdss12_m, sdss12_w = get_basic_stat_DR12(sdss_12_portSF_sa, 'Z', 'Z_ERR', 'Portsmouth SF Salpeter & SDSS & 12 ', 0.)

sys.exit()

def plotMF_raw(imf='kroupa'):
    sdss14_sel, sdss14_m, sdss14_w = get_basic_stat_FF(sdss14, 'Z', 'Z_ERR', 'ZWARNING', 'SDSS14', 0., imf=imf)
    boss14_sel, boss14_m, boss14_w = get_basic_stat_FF(boss14, 'Z_NOQSO', 'Z_ERR_NOQSO', 'ZWARNING_NOQSO', 'BOSS14', 0., imf=imf)
    sdss12_sel, sdss12_m, sdss12_w = get_basic_stat_DR12(sdss12, 'Z', 'Z_ERR', 'SDSS12', 0.)
    boss12_sel, boss12_m, boss12_w = get_basic_stat_DR12(boss12, 'Z', 'Z_ERR', 'BOSS12', 0.)

    
    dlog10m = 0.25
    mbins = n.arange(8,12.5,dlog10m)
    
    p.figure(1, (8,8))
    p.plot(mbins, smf01(10**mbins), label='Ilbert 13, 0.2<z<0.5', ls='dashed')
    
    x, y, ye = get_hist(sdss14_m[sdss14_sel], weights = sdss14_w[sdss14_sel]/(dlog10m*n.log(10)*area_sdss14*volume_per_deg2_val), mbins = mbins)
    p.errorbar(x, y, yerr = ye, label='SDSS14', lw=1)
    x, y, ye = get_hist(sdss12_m[sdss12_sel], weights = sdss12_w[sdss12_sel]/(dlog10m*n.log(10)*area_sdss12*volume_per_deg2_val), mbins = mbins)
    p.errorbar(x, y, yerr = ye, label='SDSS12', lw=1)
    x, y, ye = get_hist(boss12_m[boss12_sel], weights = boss12_w[boss12_sel]/(dlog10m*n.log(10)*area_boss12*volume_per_deg2_val), mbins = mbins)
    p.errorbar(x, y, yerr = ye, label='BOSS12', lw=0.5)
    x, y, ye = get_hist(boss14_m[boss14_sel], weights = boss14_w[boss14_sel]/(dlog10m*n.log(10)*area_boss14*volume_per_deg2_val), mbins = mbins)
    p.errorbar(x, y, yerr = ye, label='BOSS14', lw=0.5)
    p.title(str(z_min)+'<z<'+str(z_max))
    p.xlabel(r'$\log_{10}$ (stellar mass '+imf+r" / $M_\odot$ )")
    p.ylabel(r'$\Phi(M)$ [Mpc$^{-3}$ dex$^{-1}$]')
    p.yscale('log')
    p.legend(loc=0, frameon = False)
    p.ylim((1e-8, 1e-2))
    p.xlim((9.5, 12.5))
    p.grid()
    p.savefig(os.path.join(out_dir, "SDSS_SMF_"+imf+"_"+str(z_min)+'_z_'+str(z_max)+".jpg" ))
    p.clf()

import time
t0 = time.time()
id14 = n.arange(len(sdss14))
id12_2_14 = n.array([ id14[(sdss14['PLATE'] == sdss12[id12]['PLATE'])&(sdss14['MJD'] == sdss12[id12]['MJD'])&(sdss14['FIBERID'] == sdss12[id12]['FIBERID'])][0] for id12 in n.arange(5000)]) # len(sdss12)) ])
print time.time() - t0

m14_i = sdss14['stellar_mass_'+imf][id12_2_14]
m12_i = sdss12['LOGMASS'][n.arange(5000)]

ok = (m12_i >8 )&(m12_i < 13 )&(m14_i >8 )&(m14_i < 13 )
m14 = m14_i[ok]
m12 = m12_i[ok]
mms = n.arange(8,13,0.02)

outP = n.polyfit(m12, m14-m12, deg=1)
p.plot(mms, n.polyval(outP, mms), 'm')
outP = n.polyfit(m12, m14-m12, deg=2)
p.plot(mms, n.polyval(outP, mms), 'm')
p.plot(m12, m14-m12, 'b,')
p.axhline(n.mean(m14-m12), color='k')
p.xlabel('log(mass) dr12')
p.ylabel(r'$\Delta \log$(mass) dr14-dr12')
p.xlim((8,13.))
p.ylim((-2, 2))
p.grid()
p.savefig(os.path.join(out_dir, "SDSS_mass_comparison"+imf+"_"+str(z_min)+'_z_'+str(z_max)+".jpg" ))
p.show()

# for each object in the catalog, assign the corresponding stellar mass from firefly
# match in plate mjd fiberid using the SDSS DR12 as a reference i.e. make a DR12+DR14 catalog.
# fraction of objects with a more accurate stellar mass estimates : error is smaller.
# match DR14 wiht DR12
# extract mean SN per pixel in a spectrum and select objects that should have a better nass estimate.
# 

plotMF_raw(imf='salpeter')

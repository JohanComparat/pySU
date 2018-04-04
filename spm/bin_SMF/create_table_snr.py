import astropy.io.fits as fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p
import numpy as n
import os
import sys
from scipy.stats import scoreatpercentile as sc 
from scipy.interpolate import interp1d

survey = sys.argv[1]

z_min, z_max = 0., 1.6
imfs = ["Chabrier_ELODIE_", "Chabrier_MILES_", "Chabrier_STELIB_", "Kroupa_ELODIE_", "Kroupa_MILES_", "Kroupa_STELIB_",  "Salpeter_ELODIE_", "Salpeter_MILES_", "Salpeter_STELIB_" ]

z_bins = n.array([0, 0.025, 0.375, 0.7, 0.85, 1.6])

key_SNR = 'SNR_ALL'

SNR_keys = n.array([  'SNR_32_35', 'SNR_35_39', 'SNR_39_41', 'SNR_41_55', 'SNR_55_68', 'SNR_68_74', 'SNR_74_93' ])
SNR_w_min = n.array([ 32, 35, 39, 41, 55, 68, 74 ])
SNR_w_max = n.array([ 35, 39, 41, 55, 68, 74, 93 ])
wl_40 =  ((z_bins[1:]+z_bins[:-1]) * 0.5 + 1)*40.
snr_ids = n.searchsorted(SNR_w_max, wl_40)
print(SNR_keys[snr_ids])

out_dir = os.path.join(os.environ['OBS_REPO'], 'spm', 'results')

#path_2_MAG_cat = os.path.join(  os.environ['HOME'], 'SDSS', "dr14_specphot_gri.fits" )
#hd = fits.open(path_2_MAG_cat)

#path_2_sdss_cat = os.path.join(  os.environ['HOME'], 'SDSS', '26', 'catalogs', "FireFly.fits" )
#path_2_eboss_cat = os.path.join(  os.environ['HOME'], 'SDSS', 'v5_10_0', 'catalogs', "FireFly.fits" )
path_2_sdss_cat = os.path.join(  os.environ['OBS_REPO'], 'SDSS', '26', 'catalogs', "FireFly.fits" )
path_2_eboss_cat = os.path.join(  os.environ['OBS_REPO'], 'SDSS', 'v5_10_0', 'catalogs', "FireFly.fits" )

# OPENS THE CATALOGS
print("Loads catalog")
if survey =='deep2':
	deep2_dir = os.path.join(os.environ['OBS_REPO'], 'DEEP2')
	path_2_deep2_cat = os.path.join( deep2_dir, "zcat.deep2.dr4.v4.LFcatalogTC.Planck13.spm.fits" )
	catalog   = fits.open(path_2_deep2_cat)[1].data

if survey =='sdss':
	catalog   = fits.open(path_2_sdss_cat)[1].data
	z_name, z_err_name, class_name, zwarning = 'Z', 'Z_ERR', 'CLASS', 'ZWARNING'

if survey =='boss':
	catalog   = fits.open(path_2_eboss_cat)[1].data
	z_name, z_err_name, class_name, zwarning = 'Z_NOQSO', 'Z_ERR_NOQSO', 'CLASS_NOQSO', 'ZWARNING_NOQSO'


IMF = imfs[0]
prf = IMF.split('_')[0]+' & '+IMF.split('_')[1]
print(IMF, prf)
name, zflg_val, prefix = prf, 0., IMF

catalog_0 = (catalog[z_err_name] > 0.) & (catalog[z_name] > catalog[z_err_name])  & (catalog[class_name]=='GALAXY')  & (catalog[zwarning]==zflg_val) & (catalog[z_name] > z_min) & (catalog[z_name] < z_max)

catalog_zOk = catalog_0 & (catalog['SNR_ALL']>0)

converged = (catalog_zOk)&(catalog[prefix+'stellar_mass'] < 10**13. ) & (catalog[prefix+'stellar_mass'] > 10**4 )  & (catalog[prefix+'stellar_mass'] > catalog[prefix+'stellar_mass_low_1sig'] ) & (catalog[prefix+'stellar_mass'] < catalog[prefix+'stellar_mass_up_1sig'] ) 

dex04 = (converged) & (catalog[prefix+'stellar_mass'] < 10**14. ) & (catalog[prefix+'stellar_mass'] > 0 )  & (catalog[prefix+'stellar_mass'] > catalog[prefix+'stellar_mass_low_1sig'] ) & (catalog[prefix+'stellar_mass'] < catalog[prefix+'stellar_mass_up_1sig'] ) & ( - n.log10(catalog[prefix+'stellar_mass_low_1sig'])  + n.log10(catalog[prefix+'stellar_mass_up_1sig']) < 0.8 )

dex02 = (dex04) & ( - n.log10(catalog[prefix+'stellar_mass_low_1sig'])  + n.log10(catalog[prefix+'stellar_mass_up_1sig']) < 0.4 )


#target_bits

program_names = n.array(list(set( catalog['PROGRAMNAME'] )))
program_names.sort()

sourcetypes = n.array(list(set( catalog['SOURCETYPE'] )))
sourcetypes.sort()

length = lambda selection : len(selection.nonzero()[0]) 

pcs_ref = list(n.arange(0., 101, 5))
g = lambda key, s1, pcs = pcs_ref : n.hstack(( length(s1), sc(catalog[key][s1], pcs) ))

sel_pg = lambda pgr : (catalog_zOk) & (catalog['PROGRAMNAME']==pgr)

sel_st = lambda pgr : (catalog_zOk) & (catalog['SOURCETYPE']==pgr)

sel0_pg = lambda pgr : (catalog_0) & (catalog['PROGRAMNAME']==pgr)

sel0_st = lambda pgr : (catalog_0) & (catalog['SOURCETYPE']==pgr)

all_galaxies = []

tpps = []

for pg in sourcetypes:
	sel_all = sel_st(pg)
	n_all = length( sel_all )
	if n_all > 100 :
		#print(pg, n_all)
		all_galaxies.append(n_all)
		all_out = []
		for z_Min, z_Max, snr_key in zip(z_bins[:-1], z_bins[1:], SNR_keys[snr_ids]):
			s_z = sel_all &(catalog[z_name] >= z_Min) & (catalog[z_name] < z_Max)
			n_z = length(s_z)
			#print(z_Min, z_Max, n_z)
			if n_z > 0 :
				#print(n.min(catalog[snr_key][s_z]), n.max(catalog[snr_key][s_z]))
				itp = interp1d(sc(catalog[snr_key][s_z], pcs_ref), pcs_ref, kind='linear', fill_value= 100., bounds_error=False)
				#print(itp.x, itp.y)
				all_out.append( [n_z, itp(5), itp(20)] )
			else :
				all_out.append([0., -1, -1])
		all_out = n.hstack((all_out))
		tpp = pg + " & " + str(int(n_all)) + " & " + " & ".join(n.array([ str(int(el)) for el in all_out]) ) + ' \\\\ \n'
		print( tpp)
		tpps.append(tpp)
	
all_galaxies = n.array(all_galaxies)
tpps = n.array(tpps)

ids = n.argsort(all_galaxies)[::-1]

out_file = os.path.join(os.environ['OBS_REPO'], 'spm', 'results', "table_comp_"+survey+"_snr_all_sourcetype_SNR_moments.tex")
f=open(out_file, 'w')
#f.write('source type & N & \multicolumn{c}{2}{N galaxies} && \multicolumn{c}{2}{SNR ALL$>0$} & \\multicolumn{c}{2}{frefly converged} & \multicolumn{c}{2}{$\sigma_{\log_M}<0.4$} & \multicolumn{c}{2}{$\sigma_{\log_M}<0.2$} \\\\ \n')
#f.write('            &   & N & %      &             & N & % & N & % & N & %  \\\\ \n')
for jj in ids :
	f.write( tpps[jj] )

f.close()

sys.exit()


#converged = (catalog_zOk)&(catalog[prefix+'stellar_mass'] < 10**13. ) & (catalog[prefix+'stellar_mass'] > 10**4 )  & (catalog[prefix+'stellar_mass'] > catalog[prefix+'stellar_mass_low_1sig'] ) & (catalog[prefix+'stellar_mass'] < catalog[prefix+'stellar_mass_up_1sig'] ) 
#dex04 = (converged) & (catalog[prefix+'stellar_mass'] < 10**14. ) & (catalog[prefix+'stellar_mass'] > 0 )  & (catalog[prefix+'stellar_mass'] > catalog[prefix+'stellar_mass_low_1sig'] ) & (catalog[prefix+'stellar_mass'] < catalog[prefix+'stellar_mass_up_1sig'] ) & ( - n.log10(catalog[prefix+'stellar_mass_low_1sig'])  + n.log10(catalog[prefix+'stellar_mass_up_1sig']) < 0.8 )
#dex02 = (dex04) & ( - n.log10(catalog[prefix+'stellar_mass_low_1sig'])  + n.log10(catalog[prefix+'stellar_mass_up_1sig']) < 0.4 )
#m_catalog = n.log10(catalog[prefix+'stellar_mass'])
#w_catalog =  n.ones_like(catalog[prefix+'stellar_mass'])
#print(ld(catalog_zOk))
#return name + " & $"+ sld(converged)+"$ ("+str(n.round(ld(converged)/ld(catalog_zOk)*100.,1))+") & $"+ sld(dex04)+"$ ("+str(n.round(ld(dex04)/ld(catalog_zOk)*100.,1))+") & $"+ sld(dex02)+ "$ ("+str(n.round(ld(dex02)/ld(catalog_zOk)*100.,1))+r") \\\\"
##return catalog_sel, m_catalog, w_catalog




sys.exit()

for IMF in imfs :
	prf = IMF.split('_')[0]+' & '+IMF.split('_')[1]
	l2w = get_basic_stat_deep2(deep2, 'ZBEST', 'ZQUALITY', prf, 2., IMF, o2=False)
	f.write(l2w + " \n")

f.write('\\hline \n')
#l2w = get_basic_stat_DR12(boss_12_portSF_kr, 'Z', 'Z_ERR', 'Portsmouth Kroupa Star-Forming  & BOSS & 12 ', 0.) 
#f.write(l2w + " \n")
#l2w = get_basic_stat_DR12(boss_12_portPA_kr, 'Z', 'Z_ERR', 'Portsmouth Kroupa Passive  & BOSS & 12 ', 0.) 
#f.write(l2w + " \n")
#l2w = get_basic_stat_DR12(boss_12_portSF_sa, 'Z', 'Z_ERR', 'Portsmouth Salpeter Star-Forming & BOSS & 12 ', 0.) 
#f.write(l2w + " \n")
#l2w = get_basic_stat_DR12(boss_12_portPA_sa, 'Z', 'Z_ERR', 'Portsmouth Salpeter Passive & BOSS & 12 ', 0.) 
#f.write(l2w + " \n")

for IMF in imfs :
	prf = IMF.split('_')[0]+' & '+IMF.split('_')[1]
	l2w = get_basic_stat_firefly_DR14(boss, 'Z_NOQSO', 'Z_ERR_NOQSO', 'CLASS_NOQSO', 'ZWARNING_NOQSO', prf, 0., IMF)
	f.write(l2w + " \n")
	
f.write('\\hline \n')

#l2w = get_basic_stat_DR12(sdss_12_portSF_kr, 'Z', 'Z_ERR', 'Portsmouth Kroupa Star-Forming  & SDSS & 12 ', 0.) 
#f.write(l2w + " \n")
#l2w = get_basic_stat_DR12(sdss_12_portPA_kr, 'Z', 'Z_ERR', 'Portsmouth Kroupa Passive  & SDSS & 12 ', 0.) 
#f.write(l2w + " \n")
#l2w = get_basic_stat_DR12(sdss_12_portSF_sa, 'Z', 'Z_ERR', 'Portsmouth Salpeter Star-Forming & SDSS & 12 ', 0.) 
#f.write(l2w + " \n")
#l2w = get_basic_stat_DR12(sdss_12_portPA_sa, 'Z', 'Z_ERR', 'Portsmouth Salpeter Passive & SDSS & 12 ', 0.) 
#f.write(l2w + " \n")

for IMF in imfs :
	prf = IMF.split('_')[0]+' & '+IMF.split('_')[1]
	l2w = get_basic_stat_firefly_DR14(sdss, 'Z', 'Z_ERR', 'CLASS', 'ZWARNING', prf, 0., IMF)
	f.write(l2w + " \n")

f.write('\\hline \n')
f.close()
#"""
out_file = os.path.join(os.environ['OBS_REPO'], 'spm', 'results', "table_2_r.tex")
f=open(out_file, 'w')

for IMF in imfs :
	prf = IMF.split('_')[0]+' & '+IMF.split('_')[1]
	l2w = get_basic_stat_deep2(deep2, 'ZBEST', 'ZQUALITY', prf, 2., IMF, o2=True)
	f.write(l2w + " \n")

f.close()


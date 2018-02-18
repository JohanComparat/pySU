import astropy.io.fits as fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p
import numpy as n
import os
import sys
from scipy.stats import scoreatpercentile as sc 

survey = sys.argv[1]

z_min, z_max = 0., 1.6
imfs = ["Chabrier_ELODIE_", "Chabrier_MILES_", "Chabrier_STELIB_", "Kroupa_ELODIE_", "Kroupa_MILES_", "Kroupa_STELIB_",  "Salpeter_ELODIE_", "Salpeter_MILES_", "Salpeter_STELIB_" ]

out_dir = os.path.join(os.environ['OBS_REPO'], 'spm', 'results')


path_2_sdss_cat = os.path.join(  os.environ['HOME'], 'SDSS', '26', 'catalogs', "FireFly.fits" )
path_2_eboss_cat = os.path.join(  os.environ['HOME'], 'SDSS', 'v5_10_0', 'catalogs', "FireFly.fits" )

# OPENS THE CATALOGS
print("Loads catalog")
if survey =='deep2':
	deep2_dir = os.path.join(os.environ['OBS_REPO'], 'DEEP2')
	path_2_deep2_cat = os.path.join( deep2_dir, "zcat.deep2.dr4.v4.LFcatalogTC.Planck13.spm.fits" )
	deep2   = fits.open(path_2_deep2_cat)[1].data

if survey =='sdss':
	sdss   = fits.open(path_2_sdss_cat)[1].data

if survey =='sdss':
	boss   = fits.open(path_2_eboss_cat)[1].data

IMF = imfs[0]
prf = IMF.split('_')[0]+' & '+IMF.split('_')[1]
print(IMF, prf)

catalog, z_name, z_err_name, class_name, zwarning, name, zflg_val, prefix = boss, 'Z_NOQSO', 'Z_ERR_NOQSO', 'CLASS_NOQSO', 'ZWARNING_NOQSO', prf, 0., IMF

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

g = lambda key, s1, pcs = n.array([10., 25., 50., 75., 90. ]) : n.hstack(( length(s1), sc(catalog[key][s1], pcs) ))

sel_pg = lambda pgr : (catalog_zOk) & (catalog['PROGRAMNAME']==pgr)

sel_st = lambda pgr : (catalog_zOk) & (catalog['SOURCETYPE']==pgr)

sel0_pg = lambda pgr : (catalog_0) & (catalog['PROGRAMNAME']==pgr)

sel0_st = lambda pgr : (catalog_0) & (catalog['SOURCETYPE']==pgr)


out_file = os.path.join(os.environ['OBS_REPO'], 'spm', 'results', "table_comp_"+survey+"_snr_all_sourcetype_N_Nsnr_Nconv_Ndex04_Ndex02.tex")
f=open(out_file, 'w')
f.write('source type & N & N galaxies & SNR ALL$>0$ & firefly converged & err$<0.4$ & err$<0.2$ \\\\')
for pg in sourcetypes:
	f.write(pg)
	out = n.array([
	length( (catalog['SOURCETYPE']==pg)),
	length( sel0_st(pg)              ),
	length( sel_st(pg)               ), 
	length( (sel_st(pg))&(converged) ),
	length( (sel_st(pg))&(dex04)     ),
	length( (sel_st(pg))&(dex02)     )
	])
	tpp = "".join(n.array([ ' & '+str(el) for el in out]) )
	print(pg, tpp)
	f.write( tpp )
	f.write(' \\\\ \n')

f.close()


out_file = os.path.join(os.environ['OBS_REPO'], 'spm', 'results', "table_comp_"+survey+"_snr_all_programname.tex")
f=open(out_file, 'w')

for pg in program_names:
	f.write(pg)
	tpp = str( g('SNR_ALL', sel_pg(pg)) )[1:-1]
	print(pg, tpp)
	f.write( tpp )
	f.write(' \\\\ \n')

f.close()


out_file = os.path.join(os.environ['OBS_REPO'], 'spm', 'results', "table_comp_"+survey+"_snr_all_sourcetype.tex")
f=open(out_file, 'w')

for pg in sourcetypes:
	f.write(pg)
	tpp = str( g('SNR_ALL', sel_st(pg)) )[1:-1]
	print(pg, tpp)
	f.write( tpp )
	f.write(' \\\\ \n')

f.close()



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


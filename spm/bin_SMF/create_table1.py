from lib_spm import *

imfs = ["Chabrier_ELODIE_", "Chabrier_MILES_", "Chabrier_STELIB_", "Kroupa_ELODIE_", "Kroupa_MILES_", "Kroupa_STELIB_",  "Salpeter_ELODIE_", "Salpeter_MILES_", "Salpeter_STELIB_" ]
"""
out_file = os.path.join(os.environ['OBS_REPO'], 'spm', 'results', "table_1_r.tex")
f=open(out_file, 'w')

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
"""
out_file = os.path.join(os.environ['OBS_REPO'], 'spm', 'results', "table_2_r.tex")
f=open(out_file, 'w')

for IMF in imfs :
	prf = IMF.split('_')[0]+' & '+IMF.split('_')[1]
	l2w = get_basic_stat_deep2(deep2, 'ZBEST', 'ZQUALITY', prf, 2., IMF, o2=True)
	f.write(l2w + " \n")

f.close()


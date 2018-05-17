from lib_spm import *

#out_dir = os.path.join('/data42s/comparat/firefly/v1_1_0/figures', 'mass-redshift-presentation')
out_dir = os.path.join(os.environ['HOME'], 'wwwDir', 'firefly')

m_bins = n.arange(6.5, 13, 0.01)

p.figure(2, (6.5, 3.5))
p.axes([0.12,0.18,0.8,0.73])


for imf in imfs:
  stellar_mass = imf+'stellar_mass'

  redshift_reliable_boss =  (boss['CLASS_NOQSO'] == "GALAXY") & ( boss['Z_ERR_NOQSO'] > 0.0) & (boss['ZWARNING_NOQSO'] == 0) & (boss['Z_NOQSO']>0.001) & (boss['Z_NOQSO'] > boss['Z_ERR_NOQSO'] ) # (boss['SN_MEDIAN_ALL'] > 0.1 ) & 
  error_reliable_boss = (boss[stellar_mass+'_up_1sig'] > boss[stellar_mass+'_low_1sig'] ) & (boss[stellar_mass+'_up_1sig'] > 0. ) & ( boss[stellar_mass+'_low_1sig'] > 0. ) & (boss[stellar_mass+'_up_1sig'] < 1e14 ) & ( boss[stellar_mass+'_low_1sig'] < 1e14 ) 
  mass_reliable_boss_02 = (boss[stellar_mass] > 1e6 ) & ( boss[stellar_mass] < 1e14 ) & ((n.log10(boss[stellar_mass+'_up_1sig']) - n.log10(boss[stellar_mass+'_low_1sig']))/2. < 0.2 )
  mass_reliable_boss_04 = (boss[stellar_mass] > 1e6 ) & ( boss[stellar_mass] < 1e14 ) & ((n.log10(boss[stellar_mass+'_up_1sig']) - n.log10(boss[stellar_mass+'_low_1sig']))/2. < 0.4)
  ok_boss_02 = (error_reliable_boss) & (mass_reliable_boss_02) & (redshift_reliable_boss)
  ok_boss_04 = (error_reliable_boss) & (mass_reliable_boss_04) & (redshift_reliable_boss)
  print( "boss 02",len(ok_boss_02.nonzero()[0]))
  print( "boss 04",len(ok_boss_04.nonzero()[0]))
  
  #Ms_02_boss = n.log10(boss[stellar_mass][ok_boss_02])
  Ms_04_boss = n.log10(boss[stellar_mass][ok_boss_04])

  p.hist(Ms_04_boss, bins = m_bins, histtype='step', label=imf[:-1] , cumulative=True, normed=True  )

p.ylabel('normed cumulative distribution')
p.xlabel(r'\log_{10}(M/M_\odot)$')
#p.yscale('log')
p.title('eBOSS')
p.legend(loc=2, frameon = False, fontsize=11)
p.xlim((6, 12.5))
p.grid()
p.savefig(os.path.join(out_dir, "mass_distribution_eboss.png" ))
p.clf()


p.figure(2, (6.5, 3.5))
p.axes([0.12,0.18,0.8,0.73])

for imf in imfs:
  stellar_mass = imf+'stellar_mass'

  redshift_reliable_sdss =  (sdss['CLASS'] == "GALAXY")       & ( sdss['Z_ERR'] > 0.0)       & (sdss['ZWARNING'] == 0)       & (sdss['Z'] > 0.001) & (sdss['Z'] > sdss['Z_ERR'] ) # (sdss['SN_MEDIAN_ALL'] > 0.1 ) &
  error_reliable_sdss = (sdss[stellar_mass+'_up_1sig'] > sdss[stellar_mass+'_low_1sig'] ) & (sdss[stellar_mass+'_up_1sig'] > 0. ) & ( sdss[stellar_mass+'_low_1sig'] > 0. ) & (sdss[stellar_mass+'_up_1sig'] < 1e14 ) & ( sdss[stellar_mass+'_low_1sig'] < 1e14 ) 
  mass_reliable_sdss_02 = (sdss[stellar_mass] > 1e6 ) & ( sdss[stellar_mass] < 1e14 ) & ((n.log10(sdss[stellar_mass+'_up_1sig']) - n.log10(sdss[stellar_mass+'_low_1sig']))/2. < 0.2 )
  mass_reliable_sdss_04 = (sdss[stellar_mass] > 1e6 ) & ( sdss[stellar_mass] < 1e14 ) & ((n.log10(sdss[stellar_mass+'_up_1sig']) - n.log10(sdss[stellar_mass+'_low_1sig']))/2. < 0.4 )
  ok_sdss_02 = (error_reliable_sdss) & (mass_reliable_sdss_02) & (redshift_reliable_sdss)
  ok_sdss_04 = (error_reliable_sdss) & (mass_reliable_sdss_04) & (redshift_reliable_sdss)

  #Ms_02_sdss = n.log10(sdss[stellar_mass][ok_sdss_02]) 
  Ms_04_sdss = n.log10(sdss[stellar_mass][ok_sdss_04])


  p.hist(Ms_04_sdss, bins = m_bins, histtype='step', label=imf[:-1] , cumulative=True, normed=True   )

p.ylabel('normed cumulative distribution')
p.xlabel(r'\log_{10}(M/M_\odot)$')
#p.yscale('log')
p.title('SDSS')
p.legend(loc=2, frameon = False, fontsize=11)
p.xlim((6, 12.5))
p.grid()
p.savefig(os.path.join(out_dir, "mass_distribution_sdss.png" ))
p.clf()


p.figure(2, (6.5, 3.5))
p.axes([0.12,0.18,0.8,0.73])


for imf in imfs:
  stellar_mass = imf+'stellar_mass'
  z_flg = 'ZQUALITY'
  z_name = 'ZBEST'
  deep2_zOk = (deep2[z_name] > 0.001) & (deep2[z_flg]>=2.) & (deep2[z_name] < 1.7) & (deep2['SSR']>0) & (deep2['TSR']>0) & (deep2['SSR']<=1.0001) & (deep2['TSR']<=1.0001)
  deep2_sel_02 = (deep2_zOk) & (deep2[stellar_mass] < 10**14. ) & (deep2[stellar_mass] > 0. )  & (deep2[stellar_mass] >= deep2[stellar_mass+'_low_1sig'] ) & (deep2[stellar_mass] <= deep2[stellar_mass+'_up_1sig'] ) & ( - n.log10(deep2[stellar_mass+'_low_1sig'])  + n.log10(deep2[stellar_mass+'_up_1sig']) < 0.4 )
  deep2_sel_04 = (deep2_zOk) & (deep2[stellar_mass] < 10**14. ) & (deep2[stellar_mass] > 0. )  & (deep2[stellar_mass] >= deep2[stellar_mass+'_low_1sig'] ) & (deep2[stellar_mass] <= deep2[stellar_mass+'_up_1sig'] ) & ( - n.log10(deep2[stellar_mass+'_low_1sig'])  + n.log10(deep2[stellar_mass+'_up_1sig']) < 0.8 )
  Ms_02_d2 = n.log10(deep2[stellar_mass][deep2_sel_02])
  Ms_04_d2 = n.log10(deep2[stellar_mass][deep2_sel_04])
  w_deep2 = 1. / (deep2['TSR'] * deep2['SSR'])
  p.hist(Ms_04_d2, bins = m_bins, histtype='step', label=imf[:-1] , cumulative=True, normed=True   )


p.ylabel('normed cumulative distribution')
p.xlabel(r'\log_{10}(M/M_\odot)$')
#p.yscale('log')
p.title('DEEP2')
p.legend(loc=2, frameon = False, fontsize=11)
p.xlim((6, 12.5))
p.grid()
p.savefig(os.path.join(out_dir, "mass_distribution_deep2.png" ))
p.clf()






















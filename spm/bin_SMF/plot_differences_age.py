from lib_spm import *
from scipy.stats import norm
#out_dir = os.path.join('/data42s/comparat/firefly/v1_1_0/figures', 'mass-redshift-presentation')
out_dir = os.path.join(os.environ['HOME'], 'wwwDir', 'firefly')

def plotDIFF(imf_ref, imf_1, imf_2, m_bins = n.arange(-10., 10., 0.1)):
	print('eboss')
	stellar_mass = imf_ref+'stellar_mass'
	age = imf_ref+'age_lightW'
	metal = imf_ref+'metallicity_lightW'
	# selections
	redshift_reliable_boss =  (boss['CLASS_NOQSO'] == "GALAXY") & ( boss['Z_ERR_NOQSO'] > 0.0) & (boss['ZWARNING_NOQSO'] == 0) & (boss['Z_NOQSO']>0.001) & (boss['Z_NOQSO'] > boss['Z_ERR_NOQSO'] ) # (boss['SN_MEDIAN_ALL'] > 0.1 ) & 
	error_reliable_boss = (boss[stellar_mass+'_up_1sig'] > boss[stellar_mass+'_low_1sig'] ) & (boss[stellar_mass+'_up_1sig'] > 0. ) & ( boss[stellar_mass+'_low_1sig'] > 0. ) & (boss[stellar_mass+'_up_1sig'] < 1e14 ) & ( boss[stellar_mass+'_low_1sig'] < 1e14 ) 	
	mass_reliable_boss_02 = (boss[stellar_mass] > 1e6 ) & ( boss[stellar_mass] < 1e14 ) & ((n.log10(boss[stellar_mass+'_up_1sig']) - n.log10(boss[stellar_mass+'_low_1sig']))/2. < 0.2 )
	mass_reliable_boss_04 = (boss[stellar_mass] > 1e6 ) & ( boss[stellar_mass] < 1e14 ) & ((n.log10(boss[stellar_mass+'_up_1sig']) - n.log10(boss[stellar_mass+'_low_1sig']))/2. < 0.4)
	ok_boss_02 = (error_reliable_boss) & (mass_reliable_boss_02) & (redshift_reliable_boss)
	ok_boss_04 = (error_reliable_boss) & (mass_reliable_boss_04) & (redshift_reliable_boss)
	# defines quantities
	A_04_ref = boss[age][ok_boss_04]
	M_04_ref = boss[stellar_mass][ok_boss_04]
	Z_04_ref = boss[metal][ok_boss_04]
	# defines errors
	eM_04_ref = (boss[stellar_mass+'_up_1sig'][ok_boss_04]-boss[stellar_mass+'_low_1sig'][ok_boss_04])/2.
	eA_04_ref = (boss[age+'_up_1sig'][ok_boss_04]-boss[age+'_low_1sig'][ok_boss_04])/2.
	eZ_04_ref = (boss[metal+'_up_1sig'][ok_boss_04]-boss[metal+'_low_1sig'][ok_boss_04])/2.
	# quantity to compare to
	stellar_mass = imf_1+'stellar_mass'
	age = imf_1+'age_lightW'
	metal = imf_1+'metallicity_lightW'
	A_04_1 = boss[age][ok_boss_04]
	M_04_1 = boss[stellar_mass][ok_boss_04]
	Z_04_1 = boss[metal][ok_boss_04]
	eM_04_1 = (boss[stellar_mass+'_up_1sig'][ok_boss_04]-boss[stellar_mass+'_low_1sig'][ok_boss_04])/2.
	eA_04_1 = (boss[age+'_up_1sig'][ok_boss_04]-boss[age+'_low_1sig'][ok_boss_04])/2.
	eZ_04_1 = (boss[metal+'_up_1sig'][ok_boss_04]-boss[metal+'_low_1sig'][ok_boss_04])/2.
	# second comparison
	stellar_mass = imf_2+'stellar_mass'
	age = imf_2+'age_lightW'
	metal = imf_2+'metallicity_lightW'
	A_04_2 = boss[age][ok_boss_04]
	M_04_2 = boss[stellar_mass][ok_boss_04]
	Z_04_2 = boss[metal][ok_boss_04]
	eM_04_2 = (boss[stellar_mass+'_up_1sig'][ok_boss_04]-boss[stellar_mass+'_low_1sig'][ok_boss_04])/2.
	eA_04_2 = (boss[age+'_up_1sig'][ok_boss_04]-boss[age+'_low_1sig'][ok_boss_04])/2.
	eZ_04_2 = (boss[metal+'_up_1sig'][ok_boss_04]-boss[metal+'_low_1sig'][ok_boss_04])/2.
	# normalized comparison ratio
	delta_m1 = (M_04_1-M_04_ref)*(eM_04_ref**2.+eM_04_1**2.)**(-0.5)
	delta_m2 = (M_04_2-M_04_ref)*(eM_04_ref**2.+eM_04_2**2.)**(-0.5)
	# figure
	p.figure(2, (6.5, 3.5))
	p.axes([0.12,0.18,0.8,0.73])
	p.hist(delta_m1, bins=m_bins, histtype='step', label=imf_1.split('_')[1]+"-"+imf_ref.split('_')[1] , normed=True )
	p.hist(delta_m2, bins=m_bins, histtype='step', label=imf_2.split('_')[1]+"-"+imf_ref.split('_')[1] , normed=True )
	p.plot(m_bins, norm.pdf(m_bins, loc=0, scale=1), label='N(0,1)', ls='dashed')
	p.ylabel('normed distribution')
	p.xlabel(r'$(M_1-M_{ref})/\sqrt{\sigma^2_{M_1}+\sigma^2_{M_{ref}}}$')
	p.title('eBOSS '+imf_ref.split('_')[0])
	p.legend(loc=2, frameon = False, fontsize=11)
	p.xlim((-7.5, 7.5))
	p.grid()
	p.savefig(os.path.join(out_dir, "delta_M_distribution_eboss"+imf_ref.split('_')[0]+".png" ))
	p.clf()
	# normalized comparison ratio
	delta_m1 = (A_04_1-A_04_ref)*(eA_04_ref**2.+eA_04_1**2.)**(-0.5)
	delta_m2 = (A_04_2-A_04_ref)*(eA_04_ref**2.+eA_04_2**2.)**(-0.5)
	# figure
	p.figure(2, (6.5, 3.5))
	p.axes([0.12,0.18,0.8,0.73])
	bad = (delta_m1==-n.inf)|(delta_m1==n.inf)
	p.hist(delta_m1[bad==False], bins=m_bins, histtype='step', label=imf_1.split('_')[1]+"-"+imf_ref.split('_')[1] , normed=True )
	bad = (delta_m2==-n.inf)|(delta_m2==n.inf)
	p.hist(delta_m2[bad==False], bins=m_bins, histtype='step', label=imf_2.split('_')[1]+"-"+imf_ref.split('_')[1] , normed=True )
	p.plot(m_bins, norm.pdf(m_bins, loc=0, scale=1), label='N(0,1)', ls='dashed')
	p.ylabel('normed distribution')
	p.xlabel(r'$(age_1-age_{ref})/\sqrt{\sigma^2_{age_1}+\sigma^2_{age_{ref}}}$')
	#p.yscale('log')
	p.title('eBOSS '+imf_ref.split('_')[0])
	p.legend(loc=2, frameon = False, fontsize=11)
	p.xlim((-7.5, 7.5))
	p.grid()
	p.savefig(os.path.join(out_dir, "delta_A_distribution_eboss"+imf_ref.split('_')[0]+".png" ))
	p.clf()
	# normalized comparison ratio
	delta_m1 = (Z_04_1-Z_04_ref)*(eZ_04_ref**2.+eZ_04_1**2.)**(-0.5)
	delta_m2 = (Z_04_2-Z_04_ref)*(eZ_04_ref**2.+eZ_04_2**2.)**(-0.5)
	# figure
	p.figure(2, (6.5, 3.5))
	p.axes([0.12,0.18,0.8,0.73])
	bad = (delta_m1==-n.inf)|(delta_m1==n.inf)|n.isnan(delta_m1)
	p.hist(delta_m1[bad==False], bins=m_bins, histtype='step', label=imf_1.split('_')[1]+"-"+imf_ref.split('_')[1] , normed=True )
	bad = (delta_m2==-n.inf)|(delta_m2==n.inf)|n.isnan(delta_m2)
	p.hist(delta_m2[bad==False], bins=m_bins, histtype='step', label=imf_2.split('_')[1]+"-"+imf_ref.split('_')[1] , normed=True )
	p.plot(m_bins, norm.pdf(m_bins, loc=0, scale=1), label='N(0,1)', ls='dashed')
	p.ylabel('normed distribution')
	p.xlabel(r'$(Z_1-Z_{ref})/\sqrt{\sigma^2_{Z_1}+\sigma^2_{Z_{ref}}}$')
	#p.yscale('log')
	p.title('eBOSS '+imf_ref.split('_')[0])
	p.legend(loc=2, frameon = False, fontsize=11)
	p.xlim((-7.5, 7.5))
	p.grid()
	p.savefig(os.path.join(out_dir, "delta_Z_distribution_eboss"+imf_ref.split('_')[0]+".png" ))
	p.clf()

plotDIFF(	imf_ref = imfs[0], 	imf_1   = imfs[1], 	imf_2   = imfs[2], m_bins = n.arange(-10., 10., 0.1))
plotDIFF(	imf_ref = imfs[3], 	imf_1   = imfs[4], 	imf_2   = imfs[5], m_bins = n.arange(-10., 10., 0.1))
plotDIFF(	imf_ref = imfs[6], 	imf_1   = imfs[7], 	imf_2   = imfs[8], m_bins = n.arange(-10., 10., 0.1))


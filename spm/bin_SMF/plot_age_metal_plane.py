from lib_spm import *

#out_dir = os.path.join('/data42s/comparat/firefly/v1_1_0/figures', 'mass-redshift-presentation')
out_dir = os.path.join(os.environ['HOME'], 'wwwDir', 'firefly', 'az')

def plot_az(imf_ref):
	stellar_mass = imf_ref+'stellar_mass'
	age = imf_ref+'age_lightW'
	metal = imf_ref+'metallicity_lightW'

	redshift_reliable_boss =  (boss['CLASS_NOQSO'] == "GALAXY") & ( boss['Z_ERR_NOQSO'] > 0.0) & (boss['ZWARNING_NOQSO'] == 0) & (boss['Z_NOQSO']>0.001) & (boss['Z_NOQSO'] > boss['Z_ERR_NOQSO'] ) # (boss['SN_MEDIAN_ALL'] > 0.1 ) & 
	#redshift_reliable_sdss =  (sdss['CLASS'] == "GALAXY")       & ( sdss['Z_ERR'] > 0.0)       & (sdss['ZWARNING'] == 0)       & (sdss['Z'] > 0.001) & (sdss['Z'] > sdss['Z_ERR'] ) # (sdss['SN_MEDIAN_ALL'] > 0.1 ) &

	error_reliable_boss = (boss[stellar_mass+'_up_1sig'] > boss[stellar_mass+'_low_1sig'] ) & (boss[stellar_mass+'_up_1sig'] > 0. ) & ( boss[stellar_mass+'_low_1sig'] > 0. ) & (boss[stellar_mass+'_up_1sig'] < 1e14 ) & ( boss[stellar_mass+'_low_1sig'] < 1e14 ) 
	#error_reliable_sdss = (sdss[stellar_mass+'_up_1sig'] > sdss[stellar_mass+'_low_1sig'] ) & (sdss[stellar_mass+'_up_1sig'] > 0. ) & ( sdss[stellar_mass+'_low_1sig'] > 0. ) & (sdss[stellar_mass+'_up_1sig'] < 1e14 ) & ( sdss[stellar_mass+'_low_1sig'] < 1e14 ) 

	mass_reliable_boss_04 = (boss[stellar_mass] > 1e6 ) & ( boss[stellar_mass] < 1e14 ) & ((n.log10(boss[stellar_mass+'_up_1sig']) - n.log10(boss[stellar_mass+'_low_1sig']))/2. < 0.4 )
	#mass_reliable_sdss_04 = (sdss[stellar_mass] > 1e6 ) & ( sdss[stellar_mass] < 1e14 ) & ((n.log10(sdss[stellar_mass+'_up_1sig']) - n.log10(sdss[stellar_mass+'_low_1sig']))/2. < 0.4 )

	ok_boss_04 = (error_reliable_boss) & (mass_reliable_boss_04) & (redshift_reliable_boss)
	#ok_sdss_04 = (error_reliable_sdss) & (mass_reliable_sdss_04) & (redshift_reliable_sdss)

	A_04_ref = n.log10(boss[age][ok_boss_04]           )
	Z_04_ref = n.log10(boss[metal][ok_boss_04]         )

	a_bins = n.arange(6.5, 10.5, 0.2)
	z_bins = n.arange(-3,3,0.2)

	XX, YY = n.meshgrid((z_bins[1:]+z_bins[:-1])/2., 0.5*(a_bins[1:]+a_bins[:-1]))

	p.figure(0, (5.5, 4.5))
	p.axes([0.2,0.2,0.7,0.7])
	HH = n.histogram2d(Z_04_ref, A_04_ref, bins=[z_bins, a_bins])[0].T
	p.scatter(XX[HH>10], YY[HH>10], c=n.log10(HH[HH>10]), s=40, edgecolors='none', marker='s' )
	p.ylabel(r'$\log_{10}(Age/yr)$')
	p.xlabel(r'$\log_{10}(Z/Z_\odot)$')
	p.colorbar(shrink=0.7, label=r'$\log_{10}(N)$')
	p.legend(loc=0, frameon = False)
	#p.xlim((-0.05, 1.4))
	#p.ylim((6.5, 12.5))
	p.grid()
	p.title('eBOSS '+r'$\sigma_M < 0.4$')
	p.savefig(os.path.join(out_dir, "age_metallicity_"+imf_ref+"eboss_04.png" ))
	p.clf()


plot_az(imf_ref = imfs[0])
plot_az(imf_ref = imfs[1])
plot_az(imf_ref = imfs[2])
plot_az(imf_ref = imfs[3])
plot_az(imf_ref = imfs[4])
plot_az(imf_ref = imfs[5])
plot_az(imf_ref = imfs[6])
plot_az(imf_ref = imfs[7])
plot_az(imf_ref = imfs[8])

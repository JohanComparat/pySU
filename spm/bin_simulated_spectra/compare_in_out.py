"""
from scipy.integrate import quad
mass_function_chabrier = lambda logm : 0.158 * n.e**(- (logm - n.log10(0.079))**2./ (2*0.69**2.) ) * 10**logm * n.log(10)
print quad(mass_function_chabrier, -5, 5)[0] # amount of stellar mass in a cubic parsec 
mass_function_salpeter = lambda logm : 0.001 * (10**logm)**-2.35
print quad(mass_function_salpeter, -5, 5)[0] # amount of stellar mass in a cubic parsec 
"""
import os
import numpy as n
#import idlsave
#import glob
import matplotlib.pyplot as p
import astropy.io.fits as fits
import sys
spm_dir = os.path.join(os.environ['DATA_DIR'], "spm")
gama_dir = os.path.join(spm_dir, "GAMAmock")
input_dir = os.path.join(gama_dir, "inputs")
plot_dir = os.path.join(gama_dir, "plots")

out_name = os.path.join(gama_dir, "catalogs","GAMA.mock.spectra.spm.in.out.fits")
cat = fits.open(out_name)[1].data

# opens the merged catalog
#cat['input_stellar_mass']
#cat['input_age_massW']
#cat['input_metallicity_massW']

bins = n.arange(-3.75, 3.76, 0.025 )
p.figure(0, (4.5,4.5))
p.axhline(0.5, ls='dashed', c='k')
p.hist(cat['input_stellar_mass']-cat['stellar_mass_Chabrier_stelib'], label='cs',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(cat['input_stellar_mass']-cat['stellar_mass_Chabrier_miles'] , label='cm',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(cat['input_stellar_mass']-cat['stellar_mass_Chabrier_elodie'], label='ce',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(cat['input_stellar_mass']-cat['stellar_mass_Kroupa_stelib']  , label='ks',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(cat['input_stellar_mass']-cat['stellar_mass_Kroupa_miles']   , label='km',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(cat['input_stellar_mass']-cat['stellar_mass_Kroupa_elodie']  , label='ke',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(cat['input_stellar_mass']-cat['stellar_mass_Salpeter_stelib'], label='ss',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(cat['input_stellar_mass']-cat['stellar_mass_Salpeter_miles'] , label='sm',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(cat['input_stellar_mass']-cat['stellar_mass_Salpeter_elodie'], label='se',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.xlim((-2.,1.))
p.xlabel(r'$\Delta\log(M)$')
p.ylabel('normed cumulative histogram')
p.legend(loc=2, fontsize=12, frameon=False)
p.grid()
#p.show()
p.savefig(os.path.join(plot_dir, "hist_delta_log_Mass.png"))
p.clf()

bins = n.arange(-1., 1.06, 0.05 )
p.figure(0, (4.5,4.5))
p.axhline(0.5, ls='dashed', c='k')
p.hist(n.log10(cat['input_age_massW'])-n.log10(cat['age_massW_Chabrier_stelib']), label='cs',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(n.log10(cat['input_age_massW'])-n.log10(cat['age_massW_Chabrier_miles'] ), label='cm',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(n.log10(cat['input_age_massW'])-n.log10(cat['age_massW_Chabrier_elodie']), label='ce',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(n.log10(cat['input_age_massW'])-n.log10(cat['age_massW_Kroupa_stelib']  ), label='ks',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(n.log10(cat['input_age_massW'])-n.log10(cat['age_massW_Kroupa_miles']   ), label='km',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(n.log10(cat['input_age_massW'])-n.log10(cat['age_massW_Kroupa_elodie']  ), label='ke',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(n.log10(cat['input_age_massW'])-n.log10(cat['age_massW_Salpeter_stelib']), label='ss',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(n.log10(cat['input_age_massW'])-n.log10(cat['age_massW_Salpeter_miles'] ), label='sm',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(n.log10(cat['input_age_massW'])-n.log10(cat['age_massW_Salpeter_elodie']), label='se',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.xlabel(r'$\Delta\log(age)$ [mass weighted]')
p.ylabel('normed cumulative histogram')
p.legend(loc=2, fontsize=12, frameon=False)
p.grid()
#p.show()
p.savefig(os.path.join(plot_dir, "hist_delta_log_Age.png"))
p.clf()


bins = n.arange(-0.1, 0.2, 0.0005 )
p.figure(0, (4.5,4.5))
p.axhline(0.5, ls='dashed', c='k')
#p.hist(cat['spm_EBV_Chabrier_stelib'], label='cs',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(cat['spm_EBV_Chabrier_miles'] , label='cm',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
#p.hist(cat['spm_EBV_Chabrier_elodie'], label='ce',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
#p.hist(cat['spm_EBV_Kroupa_stelib']  , label='ks',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(cat['spm_EBV_Kroupa_miles']   , label='km',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
#p.hist(cat['spm_EBV_Kroupa_elodie']  , label='ke',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
#p.hist(cat['spm_EBV_Salpeter_stelib'], label='ss',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(cat['spm_EBV_Salpeter_miles'] , label='sm',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
#p.hist(cat['spm_EBV_Salpeter_elodie'], label='se',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.xlim((-0.02,0.075))
p.xlabel(r'SPM fitted $E(B-V)$')
p.ylabel('normed cumulative histogram')
p.legend(loc=4, fontsize=12, frameon=False)
p.grid()
#p.show()
p.savefig(os.path.join(plot_dir, "hist_ebv_miles.png"))
p.clf()


bins = n.arange(-0.1, 0.2, 0.0005 )
p.figure(0, (4.5,4.5))
p.axhline(0.5, ls='dashed', c='k')
#p.hist(cat['spm_EBV_Chabrier_stelib'], label='cs',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
#p.hist(cat['spm_EBV_Chabrier_miles'] , label='cm',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(cat['spm_EBV_Chabrier_elodie'], label='ce',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
#p.hist(cat['spm_EBV_Kroupa_stelib']  , label='ks',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
#p.hist(cat['spm_EBV_Kroupa_miles']   , label='km',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(cat['spm_EBV_Kroupa_elodie']  , label='ke',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
#p.hist(cat['spm_EBV_Salpeter_stelib'], label='ss',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
#p.hist(cat['spm_EBV_Salpeter_miles'] , label='sm',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(cat['spm_EBV_Salpeter_elodie'], label='se',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.xlim((-0.02,0.075))
p.xlabel(r'SPM fitted $E(B-V)$')
p.ylabel('normed cumulative histogram')
p.legend(loc=4, fontsize=12, frameon=False)
p.grid()
#p.show()
p.savefig(os.path.join(plot_dir, "hist_ebv_elodie.png"))
p.clf()

bins = n.arange(-0.1, 0.2, 0.0005 )
p.figure(0, (4.5,4.5))
p.axhline(0.5, ls='dashed', c='k')
p.hist(cat['spm_EBV_Chabrier_stelib'], label='cs',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
#p.hist(cat['spm_EBV_Chabrier_miles'] , label='cm',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
#p.hist(cat['spm_EBV_Chabrier_elodie'], label='ce',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(cat['spm_EBV_Kroupa_stelib']  , label='ks',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
#p.hist(cat['spm_EBV_Kroupa_miles']   , label='km',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
#p.hist(cat['spm_EBV_Kroupa_elodie']  , label='ke',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.hist(cat['spm_EBV_Salpeter_stelib'], label='ss',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
#p.hist(cat['spm_EBV_Salpeter_miles'] , label='sm',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
#p.hist(cat['spm_EBV_Salpeter_elodie'], label='se',rasterized=True, histtype='step', bins=bins, normed=True, cumulative = True)
p.xlim((-0.02,0.075))
p.xlabel(r'SPM fitted $E(B-V)$')
p.ylabel('normed cumulative histogram')
p.legend(loc=4, fontsize=12, frameon=False)
p.grid()
#p.show()
p.savefig(os.path.join(plot_dir, "hist_ebv_stelib.png"))
p.clf()
sys.exit()

def create_figure_add_chi2(filename, model_hdus):
	sp=gs.GalaxySpectrumFIREFLY(filename, milky_way_reddening=False)
	sp.openGAMAsimulatedSpectrum()
	spec = interp1d(sp.restframe_wavelength, sp.flux)
	err = interp1d(sp.restframe_wavelength, sp.error)
	wl_data_min = n.min(sp.restframe_wavelength)
	wl_data_max = n.max(sp.restframe_wavelength)

	# now creates the figure per model 
	fig = p.figure(0, figsize = (7, 10), frameon=False)#, tight_layout=True)
	rect = 0.2, 0.15, 0.85, 0.95
	#ax = fig.add_axes(rect, frameon=False)

	# panel with the spectrum
	fig.add_subplot(3,1,1)
	p.plot(sp.restframe_wavelength[::2], sp.flux[::2], 'k', rasterized =True, alpha=0.5, label='data')
	p.plot(model_hdus[0].data['wavelength'], model_hdus[0].data['model_flux'], label='model')
	p.yscale('log')
	mean_data = n.median(sp.flux)
	p.ylim((mean_data/8., mean_data*8.))
	p.xlabel('Wavelength [Angstrom]')
	p.ylabel(r'Flux [$f_\lambda$ $10^{-17}$ erg/cm2/s/A]')
	p.title(os.path.basename(filename))

	# second panel distribution of residuals
	fig.add_subplot(3,1,2)

	for hdu in model_hdus:
		#print hdu
		ok_model = (hdu.data['wavelength']>wl_data_min)&(hdu.data['wavelength']<wl_data_max)
		wl_model = hdu.data['wavelength'][ok_model]
		#print spec(wl_model),hdu.data['model_flux'][ok_model],err(wl_model)
		chi2s=(spec(wl_model)-hdu.data['model_flux'][ok_model])/err(wl_model)
		p.hist(chi2s, bins = n.arange(-2,2,0.1), normed = True, histtype='step', label=hdu.header['IMF']+hdu.header['library']+", EBV="+str(n.round(hdu.header['EBV'],3))+r", $\chi^2=$"+str(n.round(n.sum(chi2s**2.)/(len(chi2s)-2.),4)))
		p.ylim((-0.02,1.02))
		#p.yscale('log')
		p.xlabel('(data-model)/error')
		p.ylabel('Normed distribution')
		hdu.header['chi2'] =  n.sum(chi2s**2.) 
		hdu.header['ndof'] =  len(chi2s)-2.
		
	p.plot(n.arange(-2,2,0.005), gaussD.pdf(n.arange(-2,2,0.005)), 'k--', label=r'N(0,1)', lw=0.5)
	p.grid()
	p.legend(frameon=False, loc=0, fontsize=8)

	fig.add_subplot(3,1,3)
	tpl = n.transpose(n.array([ [
		hdu.header['age_lightW'],
		hdu.header['stellar_mass'],
		hdu.header['age_lightW_up']-hdu.header['age_lightW'], 
		hdu.header['age_lightW']-hdu.header['age_lightW_low'],
		hdu.header['stellar_mass_up']-hdu.header['stellar_mass'],
		hdu.header['stellar_mass']-hdu.header['stellar_mass_low']]
		for hdu in model_hdus ]))


	p.errorbar(tpl[0], tpl[1], xerr=[tpl[2], tpl[3]], yerr=[tpl[4], tpl[5]], barsabove=True, fmt='o')
	#p.axvline(prihdr['age_universe'], color='r', ls='dashed')
	idsUP = n.argsort(tpl[1])

	iterList = model_hdus[idsUP]
	for jj, hdu in enumerate(iterList):
		p.annotate(hdu.header['IMF']+hdu.header['library']+r", $\log(Z/Z_\odot)=$"+str(n.round(hdu.header['metallicity_lightW'],4)), 
				xy = (hdu.header['age_lightW'], hdu.header['stellar_mass']), xycoords='data', 
				xytext=(0.85, (jj+0.5)/len(iterList)), textcoords='axes fraction', 
				arrowprops=dict(facecolor='black', shrink=0.05, width=0.2, headwidth=3), 
				horizontalalignment='right', verticalalignment='top', fontsize=9)

	p.ylabel(r'$\log_{10}(M/[M_\odot])$')
	p.xlabel(r'$\log_{10}(age/[yr])$')
	#p.ylim((9,12.5))
	p.grid()
	p.savefig(path_2_im(filename))
	p.clf()
	return model_hdus

import time
t0t=time.time()
from os.path import join
import os
import numpy as n
import glob 
import sys 

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

import astropy.io.fits as fits

from scipy.interpolate import interp1d
from scipy.stats import norm as gaussD

import GalaxySpectrumFIREFLY as gs
import StellarPopulationModel as spm

def create_tbhdu(sp_cha, imf, lib):
	c1 = fits.Column(name='wavelength', format='D', unit='Angstrom', array=sp_cha[1].data['wavelength'])
	c2 = fits.Column(name='model_flux', format='D', unit='1e-17 erg/cm2/s', array=sp_cha[1].data['firefly_model'])
	coldefs = fits.ColDefs([c1, c2])
	tbhdu = fits.BinTableHDU.from_columns(coldefs)
	tbhdu.header['HIERARCH library'] = lib 
	tbhdu.header['HIERARCH IMF'] = imf 
	tbhdu.header['HIERARCH age_lightW']             = sp_cha[1].header['age_lightW_mean']                                
	tbhdu.header['HIERARCH age_lightW_up']          = sp_cha[1].header['age_lightW_mean_up']                                        
	tbhdu.header['HIERARCH age_lightW_low']         = sp_cha[1].header['age_lightW_mean_low']                                       
	tbhdu.header['HIERARCH metallicity_lightW']     = sp_cha[1].header['metallicity_lightW_mean']                                     
	tbhdu.header['HIERARCH metallicity_lightW_up']  = sp_cha[1].header['metallicity_lightW_mean_up']                                          
	tbhdu.header['HIERARCH metallicity_lightW_low'] = sp_cha[1].header['metallicity_lightW_mean_low']                                           
	tbhdu.header['HIERARCH age_massW']              = sp_cha[1].header['age_massW_mean']                                                 
	tbhdu.header['HIERARCH age_massW_up']           = sp_cha[1].header['age_massW_mean_up']                                                    
	tbhdu.header['HIERARCH age_massW_low']          = sp_cha[1].header['age_massW_mean_low']                                                   
	tbhdu.header['HIERARCH metallicity_massW']      = sp_cha[1].header['metallicity_massW_mean']                                       
	tbhdu.header['HIERARCH metallicity_massW_up']   = sp_cha[1].header['metallicity_massW_mean_up']                                           
	tbhdu.header['HIERARCH metallicity_massW_low']  = sp_cha[1].header['metallicity_massW_mean_low']                                          
	tbhdu.header['HIERARCH EBV']                    = sp_cha[1].header['EBV']                                                                  
	tbhdu.header['HIERARCH stellar_mass']           = sp_cha[1].header['stellar_mass_mean']                                                
	tbhdu.header['HIERARCH stellar_mass_up']        = sp_cha[1].header['stellar_mass_mean_up']                                                   
	tbhdu.header['HIERARCH stellar_mass_low']       = sp_cha[1].header['stellar_mass_mean_low']                                                  
	tbhdu.header['HIERARCH ssp_number']             = sp_cha[1].header['ssp_number']

	for el in sp_cha[1].header[33:]:
		tbhdu.header['HIERARCH '+el] = sp_cha[1].header[el]

	return tbhdu

env = 'DEEP2_DIR'

out_dir = os.path.join(os.environ[env], 'stellarpop')
im_dir = os.path.join(os.environ[env], 'stellarpop', 'images')

init_cat=join(os.environ['DEEP2_DIR'], "catalogs", "zcat.deep2.dr4.v4.LFcatalogTC.Planck15.fits")
summary_catalog = join(os.environ['DEEP2_DIR'], "catalogs", "zcat.deep2.dr4.v4.LFcatalogTC.Planck15.spm.fits")
hdu_orig_table = fits.open(init_cat)
orig_table = hdu_orig_table[1].data
orig_cols = orig_table.columns

catalog_entry = orig_table[0]

topdirs = join( os.environ['DEEP2_DIR'], 'stellarpop-*', 'stellarpop')

def make_summary(catalog_entry):
	tbhdus = []
	#table_all = []
	#for catalog_entry in orig_table:
	mask=str(catalog_entry['MASK'])
	objno=str(catalog_entry['OBJNO'])
	
	# defines output files 
	out_file = 'spFly-deep2-'+mask+'-'+objno+'.fits'
	path_2_out_file = os.path.join(out_dir, out_file)
	wwwwwqq	
	im_file = 'spFly-deep2-'+mask+'-'+objno+'.png'
	path_2_im_file = os.path.join(im_dir, im_file)
	
	# now gets the models
	models = n.array(glob.glob(join(topdirs, 'spFly-deep2-'+mask+'-'+objno+"*.fits")))
	models.sort()
		
	print mask, objno
	path_to_spectrum = glob.glob(join(os.environ['DEEP2_DIR'], 'spectra', mask, '*', '*' + objno + '*_fc_tc.dat'))
	if len(path_to_spectrum)>=1 and len(models)>=2 and os.path.isfile(path_2_out_file)==False:
		# open the observation file
		spe=gs.GalaxySpectrumFIREFLY("-", milky_way_reddening=True)
		spe.openObservedDEEP2pectrum(catalog_entry)
		
		spec = interp1d(spe.wavelength, spe.flux)
		err = interp1d(spe.wavelength, spe.error)
		wl_data_max = n.max(spe.wavelength)
		wl_data_min = n.min(spe.wavelength)
		N_data_points = len(spe.wavelength)

			
		for path_2_model in models:
			imf = path_2_model.split('/')[-3].split('-')[2]
			lib = path_2_model.split('/')[-3].split('-')[3]
			tb = create_tbhdu(fits.open(path_2_model), imf, lib)
			tbhdus.append(tb)

		prihdr = fits.Header()

		prihdr['file']   = out_file
		prihdr['mask']  = mask
		prihdr['objno']    = objno
		prihdr['models'] = 'Maraston_2011'
		prihdr['fitter'] = 'FIREFLY'
		#prihdr['ageMin'] = tbhdus[0].header['ageMin']
		#prihdr['ageMax'] = tbhdus[0].header['ageMax']
		#prihdr['Zmin']   = tbhdus[0].header['Zmin']
		#prihdr['Zmax']   = tbhdus[0].header['Zmax']
        #
		#prihdr['HIERARCH age_universe'] = tbhdus[1].header['age_universe']                                             
		prihdr['HIERARCH redshift']     = spe.redshift                                         


		# now creates the figure per model 
		fig = p.figure(0, figsize = (7, 10), frameon=False)#, tight_layout=True)
		rect = 0.2, 0.15, 0.85, 0.95
		#ax = fig.add_axes(rect, frameon=False)

		# panel with the spectrum
		fig.add_subplot(3,1,1)
		p.plot(spe.wavelength[::2], spe.flux[::2], 'k', rasterized =True, alpha=0.5)
		p.yscale('log')
		mean_data = n.median(spe.flux)
		p.ylim((mean_data/8., mean_data*8.))
		p.xlabel('Wavelength [Angstrom]')
		p.ylabel(r'Flux [$f_\lambda$ $10^{-17}$ erg/cm2/s/A]')
		p.title("mask=" + mask + ", objno=" + objno + ", z=" + str(n.round(spe.redshift,3)))
		# second panel distribution of residuals
		fig.add_subplot(3,1,2)

		for hdu in tbhdus:
			ok_model = (hdu.data['wavelength']>wl_data_min)&(hdu.data['wavelength']<wl_data_max)
			wl_model = hdu.data['wavelength'][ok_model]
			#p.plot(wl_model, (spec(wl_model)-hdu.data['model_flux'][ok_model])/err(wl_model), 'k', rasterized =True, alpha=0.5)
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
			for hdu in tbhdus ]))

		
		p.errorbar(tpl[0], tpl[1], xerr=[tpl[2], tpl[3]], yerr=[tpl[4], tpl[5]], barsabove=True, fmt='o')
		#p.axvline(prihdr['age_universe'], color='r', ls='dashed')
		idsUP = n.argsort(tpl[1])

		iterList = n.array(tbhdus)[idsUP]
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
		p.savefig(path_2_im_file)
		p.clf()

		prihdu = fits.PrimaryHDU(header=prihdr)
		newlist = [prihdu]
		for el in tbhdus:
			newlist.append(el)
		thdulist = fits.HDUList(newlist) # , tbhdu_cha_nd, tbhdu_kr_nd, tbhdu_sa_nd, tbhdu_cha_el, tbhdu_kr_el, tbhdu_sa_el ])
		if os.path.isfile(path_2_out_file ):
			os.remove(path_2_out_file )

		thdulist.writeto( path_2_out_file )

		print time.time()-t0t
		return 1.
	else:
		return 0.
		

for catalog_entry in orig_table[3100:]:
	print make_summary(catalog_entry)




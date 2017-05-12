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

env = 'DATA_DIR'
spec_dir = join( os.environ[env], "spm", "GAMAmock")
filenames = n.array(glob.glob(os.path.join(spec_dir, "gal_*.dat")))
filenames.sort()
print "N files=",len(filenames)

stellarpop_dir = join( os.environ[env], 
"spm", "GAMAmock", 'stellarpop')
out_dir = join( os.environ[env], "spm", "GAMAmock", 'results')
im_dir = os.path.join(os.environ[env], "spm", "GAMAmock", 'images')


suffixes = n.array(["-miles-cha.fits","-elodie-cha.fits","-telib-cha.fits","-miles-ss.fits", "-elodie-ss.fits","-stelib-ss.fits","-miles-kr.fits", "-elodie-kr.fits","-stelib-kr.fits"])
imfs = n.array(["Chabrier", "Chabrier","Chabrier","Salpeter","Salpeter","Salpeter","Kroupa","Kroupa","Kroupa"])
libs = n.array(["miles", "elodie","stelib","miles","elodie","stelib","miles","elodie","stelib"])

path_2_out = lambda filename : os.path.join(out_dir, os.path.basename(filename)[:-4]+".fits")
path_2_im = lambda filename : os.path.join(im_dir, os.path.basename(filename)[:-4]+".png")

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

def create_hdu_list(filename):
	model_hdus =  n.array([create_tbhdu(fits.open(os.path.join(stellarpop_dir, os.path.basename(filename)[:-4]+suff)), imf, lib) for imf, lib, suff in zip(imfs, libs, suffixes)])
	print "N hdus=",len(model_hdus)
	#print model_hdus[0].header				 
	prihdr = fits.Header()
	prihdr['file']   = os.path.basename(filename)[:-4]
	prihdr['models'] = 'Maraston_2011'
	prihdr['fitter'] = 'FIREFLY'
	return prihdr, model_hdus

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

def write_summary(filename):
	t0t = time.time()
	prihdr, model_hdus = create_hdu_list(filename)
	bla = create_figure_add_chi2(filename, model_hdus)
	model_hdus_list = bla.tolist()
	prihdu = fits.PrimaryHDU(header=prihdr)
	model_hdus_list.insert(0,prihdu)
	thdulist = fits.HDUList( model_hdus_list )
	path_2_out_file = path_2_out(filename)
	if os.path.isfile(path_2_out_file ):
		os.remove(path_2_out_file )

	thdulist.writeto( path_2_out_file )

	print time.time()-t0t, "seconds"

for filename in filenames[::-1]:
	write_summary(filename)




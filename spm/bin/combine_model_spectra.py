import time
t0t=time.time()
from os.path import join
import os
import numpy as n
import glob 
import sys 

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p

import astropy.io.fits as fits

from scipy.interpolate import interp1d
from scipy.stats import norm as gaussD

plate   = sys.argv[1]
mjd     = sys.argv[2] 
fiberid = sys.argv[3] 
env = sys.argv[4]
#env = 'EBOSSDR14_DIR'

if env == 'EBOSSDR14_DIR':
	z_name = 'Z_NOQSO' 
if env == 'SDSSDR12_DIR':
	z_name = 'Z' 

# open the observation file
obs = fits.open(os.path.join(os.environ[env], 'spectra', plate, 'spec-'+plate+'-'+mjd+'-'+fiberid+".fits"))
wl_data = 10**obs[1].data['loglam']/(1+obs[2].data[z_name])
fl_data = obs[1].data['flux']
err_data = obs[1].data['ivar']**(-0.5)
ok_data = (obs[1].data['ivar']>0)
spec = interp1d(wl_data[ok_data], fl_data[ok_data])
err = interp1d(wl_data[ok_data], err_data[ok_data])
wl_data_max = n.max(wl_data[ok_data])
wl_data_min = n.min(wl_data[ok_data])
N_data_points = len(wl_data)


#dirs = ['stellarpop-m11-salpeter', 'stellarpop-m11-kroupa', 'stellarpop-m11-chabrier', 'stellarpop-m11-salpeter-stelib', 'stellarpop-m11-kroupa-stelib', 'stellarpop-m11-chabrier-stelib', 'stellarpop-m11-salpeter-elodie', 'stellarpop-m11-kroupa-elodie', 'stellarpop-m11-chabrier-elodie'] 
#suffixs = ["-ss.fits", "-kr.fits", "-cha.fits", "-ss.fits", "-kr.fits", "-cha.fits", "-ss.fits", "-kr.fits", "-cha.fits"]

dirs = ['stellarpop-m11-chabrier', 'stellarpop-m11-kroupa', 'stellarpop-m11-salpeter'] 
suffixs = ["-cha.fits", "-kr.fits", "-ss.fits"]

print plate, mjd, fiberid

sp_cha = fits.open(os.path.join(os.environ[env], dirs[0], 'stellarpop', plate, 'spFly-'+plate+'-'+mjd+'-'+fiberid+suffixs[0]))
sp_kr  = fits.open(os.path.join(os.environ[env], dirs[1], 'stellarpop', plate, 'spFly-'+plate+'-'+mjd+'-'+fiberid+suffixs[1]))
sp_sa  = fits.open(os.path.join(os.environ[env], dirs[2], 'stellarpop', plate, 'spFly-'+plate+'-'+mjd+'-'+fiberid+suffixs[2]))
#sp_cha_nd = fits.open(os.path.join(os.environ[env],dirs[3], 'stellarpop', plate, 'spFly-'+plate+'-'+mjd+'-'+fiberid+suffixs[3]))
#sp_kr_nd = fits.open(os.path.join(os.environ[env], dirs[4], 'stellarpop', plate, 'spFly-'+plate+'-'+mjd+'-'+fiberid+suffixs[4]))
#sp_sa_nd = fits.open(os.path.join(os.environ[env], dirs[5], 'stellarpop', plate, 'spFly-'+plate+'-'+mjd+'-'+fiberid+suffixs[5]))
#sp_cha_el = fits.open(os.path.join(os.environ[env],dirs[6], 'stellarpop', plate, 'spFly-'+plate+'-'+mjd+'-'+fiberid+suffixs[6]))
#sp_kr_el = fits.open(os.path.join(os.environ[env], dirs[7], 'stellarpop', plate, 'spFly-'+plate+'-'+mjd+'-'+fiberid+suffixs[7]))
#sp_sa_el = fits.open(os.path.join(os.environ[env], dirs[8], 'stellarpop', plate, 'spFly-'+plate+'-'+mjd+'-'+fiberid+suffixs[8]))

out_dir = os.path.join(os.environ[env], 'stellarpop', plate)
im_dir = os.path.join(os.environ[env], 'stellarpop', plate, 'images')

if os.path.isdir(out_dir)==False:
	os.makedirs(out_dir)
if os.path.isdir(im_dir)==False:
	os.makedirs(im_dir)
	
out_file = 'spFly-'+plate+'-'+mjd+'-'+fiberid+'.fits'
path_2_out_file = os.path.join(out_dir, out_file)

im_file = 'spFly-'+plate+'-'+mjd+'-'+fiberid+'.png'
path_2_im_file = os.path.join(im_dir, im_file)


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

tbhdu_cha = create_tbhdu(sp_cha, 'Chabrier', 'MILES')
tbhdu_kr  = create_tbhdu(sp_kr, 'Kroupa'   , 'MILES')
tbhdu_sa  = create_tbhdu(sp_sa, 'Salpeter' , 'MILES')
#tbhdu_cha_nd = create_tbhdu(sp_cha_nd, 'Chabrier', 'STELIB')
#tbhdu_kr_nd  = create_tbhdu(sp_kr_nd, 'Kroupa'   , 'STELIB')
#tbhdu_sa_nd  = create_tbhdu(sp_sa_nd, 'Salpeter' , 'STELIB')
#tbhdu_cha_el = create_tbhdu(sp_cha_el, 'Chabrier', 'ELODIE')
#tbhdu_kr_el  = create_tbhdu(sp_kr_el, 'Kroupa'   , 'ELODIE')
#tbhdu_sa_el  = create_tbhdu(sp_sa_el, 'Salpeter' , 'ELODIE')


prihdr = fits.Header()

prihdr['file']   = out_file
prihdr['plate']  = int(plate)
prihdr['mjd']    = int(mjd)
prihdr['fiberid']= int(fiberid)
prihdr['models'] = 'Maraston_2011'
prihdr['fitter'] = 'FIREFLY'
prihdr['model']  = sp_cha[0].header['model']
prihdr['ageMin'] = sp_cha[0].header['ageMin']
prihdr['ageMax'] = sp_cha[0].header['ageMax']
prihdr['Zmin']   = sp_cha[0].header['Zmin']
prihdr['Zmax']   = sp_cha[0].header['Zmax']

prihdr['HIERARCH age_universe'] = sp_cha[1].header['HIERARCH age_universe']                                             
prihdr['HIERARCH redshift']     = sp_cha[1].header['HIERARCH redshift']                                         


# now creates the figure per model 
fig = p.figure(0, figsize = (7, 10), frameon=False)#, tight_layout=True)
rect = 0.2, 0.15, 0.85, 0.95
#ax = fig.add_axes(rect, frameon=False)

# panel with the spectrum
fig.add_subplot(3,1,1)
p.plot(wl_data[::2], fl_data[::2], 'k', rasterized =True, alpha=0.5)
p.yscale('log')
mean_data = n.median(fl_data)
p.ylim((mean_data/8., mean_data*8.))
p.xlabel('Wavelength [Angstrom]')
p.ylabel(r'Flux [$f_\lambda$ $10^{-17}$ erg/cm2/s/A]')
p.title("plate=" + plate + ", mjd=" + mjd + ", fiber=" + fiberid + ", z=" + str(n.round(obs[2].data[z_name][0],3)))
# second panel distribution of residuals
fig.add_subplot(3,1,2)

for hdu in [tbhdu_cha, tbhdu_kr, tbhdu_sa]:
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
	for hdu in [tbhdu_cha, tbhdu_kr, tbhdu_sa] ]))

 
p.errorbar(tpl[0], tpl[1], xerr=[tpl[2], tpl[3]], yerr=[tpl[4], tpl[5]], barsabove=True, fmt='o')
#p.axvline(prihdr['age_universe'], color='r', ls='dashed')
idsUP = n.argsort(tpl[1])

iterList = n.array([tbhdu_cha, tbhdu_kr, tbhdu_sa])[idsUP]
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

thdulist = fits.HDUList([prihdu, tbhdu_cha, tbhdu_kr, tbhdu_sa]) # , tbhdu_cha_nd, tbhdu_kr_nd, tbhdu_sa_nd, tbhdu_cha_el, tbhdu_kr_el, tbhdu_sa_el ])
if os.path.isfile(path_2_out_file ):
	os.remove(path_2_out_file )

thdulist.writeto( path_2_out_file )

print time.time()-t0t






import time
t0t=time.time()
from os.path import join
import os
import numpy as n
import glob 
import sys 
import astropy.io.fits as fits
import matplotlib.pyplot as p
from scipy.interpolate import interp1d
from scipy.stats import norm as gaussD
plate   = sys.argv[1]
mjd     = sys.argv[2] 
fiberid = sys.argv[3] 

env = 'EBOSSDR14_DIR'

dirs = ['stellarpop'] 
suffixs = [".fits"]

print plate, mjd, fiberid

obs = fits.open(os.path.join(os.environ[env], 'spectra', plate, 'spec-'+plate+'-'+mjd+'-'+fiberid+suffixs[0]))
wl_data = 10**obs[1].data['loglam']/(1+obs[2].data['Z'])
fl_data = obs[1].data['flux']
err_data = obs[1].data['ivar']**(-0.5)
ok_data = (obs[1].data['ivar']>0)
spec = interp1d(wl_data[ok_data], fl_data[ok_data])
err = interp1d(wl_data[ok_data], err_data[ok_data])
wl_data_max = n.max(wl_data[ok_data])
wl_data_min = n.min(wl_data[ok_data])
N_data_points = len(wl_data)

hdus = fits.open(os.path.join(os.environ[env], dirs[0], plate, 'spFly-'+plate+'-'+mjd+'-'+fiberid+suffixs[0]))

out_dir = os.path.join(os.environ[env], 'stellarpop_plots', plate)
if os.path.isdir(out_dir)==False:
	os.makedirs(out_dir)

out_file = 'flyPlot-'+plate+'-'+mjd+'-'+fiberid+'.png'
path_2_out_file = os.path.join(out_dir, out_file)


#def plot_spm(hdus):
    
Npanels = 3 # len(hdus)

fig = p.figure(0, figsize = (7, 10), frameon=False)#, tight_layout=True)
rect = 0.2, 0.15, 0.85, 0.95
#ax = fig.add_axes(rect, frameon=False)

# panel with the spectrum
fig.add_subplot(Npanels,1,1)
p.plot(wl_data[::2], fl_data[::2], 'k', rasterized =True, alpha=0.5)
p.yscale('log')
mean_data = n.median(fl_data)
p.ylim((mean_data/8., mean_data*8.))
p.xlabel('Wavelength [Angstrom]')
p.ylabel(r'Flux [$f_\lambda$ $10^{-17}$ erg/cm2/s/A]')
p.title("plate=" + plate + ", mjd=" + mjd + ", fiber=" + fiberid + ", z=" + str(n.round(obs[2].data['Z'][0],4)))
# second panel distribution of residuals
fig.add_subplot(Npanels,1,2)
	
for ii in range(1, len(hdus)):
	ok_model = (hdus[ii].data['wavelength']>wl_data_min)&(hdus[ii].data['wavelength']<wl_data_max)
	wl_model = hdus[ii].data['wavelength'][ok_model]
	#p.plot(wl_model, (spec(wl_model)-hdus[ii].data['model_flux'][ok_model])/err(wl_model), 'k', rasterized =True, alpha=0.5)
	chi2s=(spec(wl_model)-hdus[ii].data['model_flux'][ok_model])/err(wl_model)
	p.hist(chi2s, bins = n.arange(-2,2,0.1), normed = True, histtype='step', label=hdus[ii].header['IMF']+hdus[ii].header['library']+", EBV="+str(n.round(hdus[ii].header['EBV'],3))+r", $\chi^2=$"+str(n.round(n.sum(chi2s**2.)/(len(chi2s)-2.),4)))
	p.ylim((-0.02,1.02))
	#p.yscale('log')
	p.xlabel('(data-model)/error')
	p.ylabel('Normed distribution')

p.plot(n.arange(-2,2,0.005), gaussD.pdf(n.arange(-2,2,0.005)), 'k--', label=r'N(0,1)', lw=0.5)
p.grid()
p.legend(frameon=False, loc=0, fontsize=8)

fig.add_subplot(Npanels,1,3)
tpl = n.transpose(n.array([ [
	hdus[ii].header['age_lightW'],
	hdus[ii].header['stellar_mass'],
	hdus[ii].header['age_lightW_up']-hdus[ii].header['age_lightW'], 
	hdus[ii].header['age_lightW']-hdus[ii].header['age_lightW_low'],
	hdus[ii].header['stellar_mass_up']-hdus[ii].header['stellar_mass'],
	hdus[ii].header['stellar_mass']-hdus[ii].header['stellar_mass_low']]
	for ii in range(1, len(hdus)) ]))

 
p.errorbar(tpl[0], tpl[1], xerr=[tpl[2], tpl[3]], yerr=[tpl[4], tpl[5]], barsabove=True, fmt='o')
p.axvline(hdus[0].header['age_universe'], color='r', ls='dashed')
idsUP = n.argsort(tpl[1])
for jj, ii in enumerate(n.arange(1, len(hdus))[idsUP]):
	p.annotate(hdus[ii].header['IMF']+hdus[ii].header['library']+", EBV="+str(n.round(hdus[ii].header['EBV'],3)), xy = (tpl[0][ii-1], tpl[1][ii-1]), xycoords='data', xytext=(0.85, (jj+0.5)*1./len(hdus)), textcoords='axes fraction', arrowprops=dict(facecolor='black', shrink=0.05, width=0.2, headwidth=3), horizontalalignment='right', verticalalignment='top')

p.ylabel(r'$log_{10}(M/M_\odot)$')
p.xlabel(r'$log_{10}(age/yr)$')
#p.ylim((9,12.5))
p.grid()
p.savefig(path_2_out_file)
p.clf()
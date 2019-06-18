import astropy.io.fits as fits
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p
import numpy as n

from cycler import cycler
# 1. Setting prop cycle on default rc parameter
p.rc('lines', linewidth=0.7)
p.rc('axes', prop_cycle=cycler('color', ["#638bd9", "#e586b6", "#dd7c25"]) )#, "#598664", "#9631a9", "#cff159", "#534f55", "#ab1519", "#89dbde"]) )

import os, glob, sys

from scipy.stats import norm

m_bins = n.arange(-4,4,0.1)

stack_dir = os.path.join(os.environ['HOME'], 'data', 'SDSS', 'eBOSS', 'ELG', 'stacks')
ff_dir = os.path.join(os.environ['HOME'], 'data', 'SDSS', 'eBOSS', 'ELG', 'stacks', 'firefly')

l_dict = {}
l_dict["O2"] = r"$[O_\mathrm{II}] $"
l_dict["O3"] = r"$[O_\mathrm{III}]$"
l_dict["Hb"] = r"$[H_\mathrm{\beta}]$"

line = sys.argv[1] # "O2"

stack_list = n.array(glob.glob(os.path.join(stack_dir, '*'+line+'*.stack')))
ascii_list = n.array(glob.glob(os.path.join(stack_dir, '*'+line+'*.ascii')))
ff_list = n.array(glob.glob(os.path.join(ff_dir, '*'+line+'*.stack')))
stack_list.sort()
ff_list.sort()

N_spec = len(stack_list)

fig_dir = os.path.join(os.environ['GIT_PYSU'], 'figures/sdss/elg/stacks/LF/',line)
if os.path.isdir(fig_dir)==False:
	os.system('mkdir '+fig_dir)

def get_data(path_2_spec, path_2_ff):
	s = fits.open(path_2_spec)
	d = fits.open(path_2_ff)
	selection_stack = ( s[1].data['NspectraPerPixel'] > n.max(s[1].data['NspectraPerPixel'])*0.5 ) & ( s[1].data['medianStack'] > 0 )
	x_data = s[1].data['wavelength']        [selection_stack]
	y_data = s[1].data['medianStack']       [selection_stack]
	y_err = s[1].data['jackknifStackErrors'][selection_stack]
	selection_model = (d[2].data['wavelength'] >2300 ) & (d[2].data['wavelength']<4000)
	x_model = d[2].data['wavelength']    [selection_model]
	y_model = d[2].data['firefly_model'] [selection_model]
	return x_data, y_data, y_err, x_model, y_model, d


for jj, (path_2_spec, path_2_ff, ascii_file) in enumerate( zip (stack_list, ff_list, ascii_list) ):
	DATA = n.loadtxt(ascii_file, unpack=True)
	z_mean = n.mean(DATA[3])
	bn_spec = os.path.basename(path_2_spec)
	bn_spec_split = n.hstack(( n.array([ el.split('_') for el in  bn_spec[:-6].split('-') ]) ))  
	path_2_figure = os.path.join( fig_dir, bn_spec+'.png')
	# A4 figure
	fig = p.figure(0, (10., 6.), frameon=False )
	print(path_2_spec, path_2_ff)
	x_data, y_data, y_err, x_model, y_model, d = get_data(path_2_spec, path_2_ff)
	# panel top left: spectrum + models
	#if jj < N_spec-1 :
	ax1=fig.add_subplot(1,1,1)  
	#ax1.set_title(bn_spec) 
	ax1.set_xlabel( 'wavelength [Angstrom, rest frame]' ) 
	ax1.set_ylabel( r'$f_\lambda [10^{-17}$ erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
	ax1.set_ylim ( (( n.min(y_model)*0.9,1.3* n.max(y_model) )) )
	ax1.set_xlim ( (2300, 4000 ) )
	#ax1.set_yscale ( 'log' )

	#if jj == N_spec -1 :
		#fig.add_subplot(N_spec,1,jj) 
		#ax1=fig.add_subplot(N_spec,1,jj+1)  
		#ax1.set_title(os.path.basename(path_2_spec)) 
		#ax1.set_xlabel( 'wavelength [Angstrom, rest frame]' ) 
		#ax1.set_ylabel( r'$f_\lambda [10^{-17}$ erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
		#ax1.set_ylim ( (( n.min(y_data)*0.9,1.1* n.max(y_data) )) )
		#ax1.set_xlim ( (( n.min(x_data)*0.9,1.1* n.max(x_data) )) )
	ii=2
	hduSPM = d[ii]
	age_log10 =  str(n.round(n.log10( 1e9*10**hduSPM.header['age_massW'] ),2))                 
	z_log10 =  str(n.round(10**hduSPM.header['metallicity_massW'],4))                 
	m_log10 =  str(n.round(hduSPM.header['stellar_mass'],2))                 
	label = d[ii].header['IMF']+" "+d[ii].header['MODEL']+r", N$_{SSP}$="+str(d[ii].header['ssp_number'])+r', $\log_{10}(age/yr)=$'+age_log10 +r', $Z/Z_\odot=$'+z_log10 +r', $\log_{10}(M/M_\odot)=$'+m_log10
	p.plot(d[ii].data['wavelength'], d[ii].data['firefly_model'], label=label)

	label_data = 'stack, N='+bn_spec_split[8]+', '+bn_spec_split[1]+r'$<\log_{10}(L$'+l_dict[line]+'/[erg/s])<'+bn_spec_split[3]+', '+bn_spec_split[4]+'<z<'+bn_spec_split[6]+r', $\bar{z}=$'+str(n.round(z_mean,2)) 
	p.plot(x_data, y_data, color='grey', lw=1., label=label_data)

	p.legend(loc=2, fontsize=9)#, frameon=False)
	p.grid()

	p.tight_layout()
	p.savefig(path_2_figure)
	p.clf()


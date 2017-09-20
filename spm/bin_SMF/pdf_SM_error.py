"""
Plots stellar mass error and stellar mass vs. mass and redshift

"""
import sys
import os
import numpy as n
import astropy.io.fits as fits
import matplotlib
matplotlib.rcParams['agg.path.chunksize'] = 2000000
matplotlib.rcParams.update({'font.size': 13})
#matplotlib.use('Agg')
import matplotlib.pyplot as p
from cycler import cycler
from lib_spm import *

out_dir = os.path.join(os.environ['OBS_REPO'], 'spm', 'results')

imfs = ["Chabrier_ELODIE_", "Chabrier_MILES_", "Chabrier_STELIB_", "Kroupa_ELODIE_", "Kroupa_MILES_", "Kroupa_STELIB_",  "Salpeter_ELODIE_", "Salpeter_MILES_", "Salpeter_STELIB_" ]

z_bins = n.arange(0, 4.1, 0.1)
m_bins = n.arange(0,14,0.5)
sn_bins = n.array([0, 0.5, 1, 2, 10, 100]) # n.logspace(-1.,2,6)
err_bins = n.hstack(( n.arange(0., 1., 0.15), n.arange(1., 5., 0.3), n.array([5., 10., 10000.]) )) 
x_err_bins = (err_bins[1:] + err_bins[:-1])/2.

prefix = 'BOSS'
imf = imfs[0] #'Chabrier'
stellar_mass = imf+'stellar_mass'

redshift_reliable =  (boss['CLASS_NOQSO'] == "GALAXY") & (boss['Z_NOQSO'] >= 0) & ( boss['Z_ERR_NOQSO'] >= 0) & (boss['ZWARNING_NOQSO'] == 0) & (boss['Z_NOQSO'] > boss['Z_ERR_NOQSO'] ) # (boss['SN_MEDIAN_ALL'] > 0.1 ) & 

error_reliable = (boss[stellar_mass+'_up'] > boss[stellar_mass+'_low'] ) & (boss[stellar_mass+'_up'] > 0. ) & ( boss[stellar_mass+'_low'] > 0. ) & (boss[stellar_mass+'_up'] < 10. ) & ( boss[stellar_mass+'_low'] < 10. ) 

mass_reliable = (boss[stellar_mass] > 0 ) & ( boss[stellar_mass] < 14. ) # & ( abs(boss[stellar_mass + '_up'] - boss[stellar_mass]) < 0.4 ) & ( abs(boss[stellar_mass + '_low'] - boss[stellar_mass]) < 0.4 )

ok = (error_reliable) & (mass_reliable) & (redshift_reliable)

snr = boss['SN_MEDIAN_ALL'][ok]
zz = boss['Z_NOQSO'][ok]
Ms = boss[stellar_mass][ok]
err_moyenne = abs((boss[stellar_mass + '_up'][ok] - boss[stellar_mass + '_low'][ok])/2.)

def plot_errPDF(degree):
	p.figure(1, (4.5, 4.5))
	p.axes([0.2,0.2,0.7,0.7])
	cc=cycler('color', ['r', 'g', 'b', 'm', 'k'])
	print "DEGREE", degree
	for ii, (sn, c) in enumerate(zip(sn_bins[:-1],cc)):
		out, xxx = n.histogram(err_moyenne[(snr>sn_bins[ii])&(snr<sn_bins[ii+1])], bins = err_bins, normed=True)
		outN = n.histogram(err_moyenne[(snr>sn_bins[ii])&(snr<sn_bins[ii+1])], bins = err_bins)[0]
		N100 = (outN>10)
		if len(out[N100])>2:
			print sn, "-", sn_bins[ii+1], " & " ,
			eb = p.errorbar(x_err_bins[N100], out[N100], xerr=(xxx[1:][N100]-xxx[:-1][N100])/2., yerr = out[N100]*outN[N100]**(-0.5), label=str(sn)+"-"+str(sn_bins[ii+1]), fmt='+', color=c['color'])
			
			outP = n.polyfit(x_err_bins[N100], n.log10(out[N100]), deg=degree, w = 1/( outN[N100]**(-0.5)) )
			x_poly =n.arange(0,3,0.05)
			p.plot(x_poly, 10**n.polyval(outP, x_poly), c=c['color'], ls='dashed', lw=0.5)
			print n.round(outP,3) , " \\\\"

	p.legend(loc=0, frameon=False)
	p.yscale('log')
	p.ylabel('pdf')
	p.xlim((-0.01, 1.5))
	p.ylim((1e-3, 10.))
	p.title("degree="+str(degree))
	p.grid()
	p.xlabel(r'$\Delta \log M$ ('+imf+')')
	p.savefig(os.path.join(out_dir, "pdf_DELTA_logM_SNR_"+imf+"_"+str(degree)+"_"+".jpg" ))
	p.clf()


#plot_errPDF(1)
#plot_errPDF(2)
#plot_errPDF(3)
plot_errPDF(4)
plot_errPDF(5)
sys.exit()



def plot_all_prognames(hdus=hdus, imf=imf, prefix=prefix, out_dir = out_dir, redshift_reliable=redshift_reliable ) :
	stellar_mass = imf+'_stellar_mass'
	out_dir = os.path.join(os.environ['DATA_DIR'], 'spm', 'results', 'catalogs', imf)

	mass_reliable = (boss[stellar_mass] > 0 ) & ( boss[stellar_mass] < 13. ) & ( abs(boss[stellar_mass + '_up'] - boss[stellar_mass]) < 0.4 ) & ( abs(boss[stellar_mass + '_low'] - boss[stellar_mass]) < 0.4 )

	#good_plates = (boss['PLATEQUALITY']=='good') &(boss['TARGETTYPE']=='science')

	all_names = set(boss['PROGRAMNAME'])
	all_names_arr = n.array(list(all_names))

	for ii in range(len(all_names_arr)):
		selection = (mass_reliable) & (redshift_reliable) & (boss['PROGRAMNAME']==all_names_arr[ii])
		N_occ = len(selection.nonzero()[0])
		print all_names_arr[ii], N_occ
		if N_occ>1:
			p.figure(1, (4.5, 4.5))
			p.axes([0.2,0.2,0.7,0.7])
			p.plot(n.log10(1.+boss['Z'][selection]), boss[stellar_mass][selection], 'k+', rasterized=True, alpha=0.5) #, label=all_names_arr[ii]
			p.ylabel(r'$\log_{10}$ (stellar mass '+imf+r" / $M_\odot$ )")
			p.axvline(n.log10(3.), ls='dashed', label='z=0.1, 0.5, 1, 2')
			p.axvline(n.log10(2.), ls='dashed')#, label='z=1')
			p.axvline(n.log10(1.1),ls='dashed')#, label='z=0.1')
			p.axvline(n.log10(1.5),ls='dashed')#, label='z=0.5')
			p.xlabel(r'$\log_{10}(1+z)$')
			#p.xscale('log')
			p.legend(loc=0, frameon = False)
			p.xlim((0.0, 0.7))
			p.ylim((6.5, 12.5))
			p.grid()
			p.title('N='+str(N_occ))
			p.savefig(os.path.join(out_dir, prefix+"_"+all_names_arr[ii]+"_redshift_mass_"+imf+".jpg" ))
			p.clf()

			p.figure(2, (4.5, 4.5))
			p.axes([0.2,0.2,0.7,0.7])
			#p.subplot(111, projection="mollweide")
			#p.plot((boss['PLUG_RA'][selection]-180.)*n.pi/180., boss['PLUG_DEC'][selection]*n.pi/180., 'k+', rasterized=True) # , label=all_names_arr[ii]
			p.plot(boss['PLUG_RA'][selection], boss['PLUG_DEC'][selection], 'k+', rasterized=True) # , label=all_names_arr[ii]
			p.title(all_names_arr[ii])#+', Ngal='+str(N_occ))
			p.xlim((0.0, 360.))
			p.ylim((-20., 85.))
			p.xlabel(r'ra [deg]')
			p.ylabel(r'dec [deg]')
			#p.legend(loc=0, frameon = False)
			p.grid()
			p.savefig(os.path.join(out_dir, prefix+"_"+all_names_arr[ii]+"_ra_dec_"+imf+".jpg" ))
			p.clf()



plot_all_prognames(hdus=hdus, imf=imf, prefix=prefix, redshift_reliable=redshift_reliable )
"""
imf = 'Salpeter'
plot_all_prognames(hdus=hdus, imf=imf, prefix=prefix,  redshift_reliable=redshift_reliable )

imf = 'Kroupa'
plot_all_prognames(hdus=hdus, imf=imf, prefix=prefix, redshift_reliable=redshift_reliable )
"""
prefix = 'SDSS'
redshift_reliable = (sdss['SN_MEDIAN_ALL'] > 0.1 ) & (sdss['CLASS'] == "GALAXY") & (sdss['Z'] >= 0) & ( sdss['Z_ERR'] >= 0) & (sdss['ZWARNING'] == 0) & (sdss['Z'] > sdss['Z_ERR'] )

imf = 'Chabrier'
plot_all_prognames(hdus=hdus, imf=imf, prefix=prefix, redshift_reliable=redshift_reliable )
"""
imf = 'Salpeter'
plot_all_prognames(hdus=hdus, imf=imf, prefix=prefix, redshift_reliable=redshift_reliable )

imf = 'Kroupa'
plot_all_prognames(hdus=hdus, imf=imf, prefix=prefix, redshift_reliable=redshift_reliable )
"""




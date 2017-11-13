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
import cPickle

out_dir = os.path.join(os.environ['OBS_REPO'], 'spm', 'results', 'convergence')

imfs = ["Chabrier_ELODIE_", "Chabrier_MILES_", "Chabrier_STELIB_", "Kroupa_ELODIE_", "Kroupa_MILES_", "Kroupa_STELIB_",  "Salpeter_ELODIE_", "Salpeter_MILES_", "Salpeter_STELIB_" ]

z_bins = n.arange(0, 4.1, 0.1)
m_bins = n.arange(0,14,0.5)
sn_bins = n.array([0, 0.5, 1, 2, 10, 100]) # n.logspace(-1.,2,6)
err_bins = n.hstack(( n.arange(0., 1., 0.15), n.arange(1., 5., 0.3), n.array([5., 10., 10000.]) )) 
x_err_bins = (err_bins[1:] + err_bins[:-1])/2.
x_err_width = (err_bins[1:] - err_bins[:-1])/2.

def get_err_distribution(boss, sdss, imf = imfs[0]):
	"""
	:param boss: data table
	:param imf: specifies the prefix to get the stellar mass in the table, to be taken from the imfs array above
	"""
	stellar_mass = imf+'stellar_mass'
	
	redshift_reliable_boss =  (boss['CLASS_NOQSO'] == "GALAXY") & ( boss['Z_ERR_NOQSO'] > 0) & (boss['ZWARNING_NOQSO'] == 0) & (boss['Z_NOQSO'] > boss['Z_ERR_NOQSO'] ) # (boss['SN_MEDIAN_ALL'] > 0.1 ) & 
	redshift_reliable_sdss =  (sdss['CLASS'] == "GALAXY")       & ( sdss['Z_ERR'] > 0)       & (sdss['ZWARNING'] == 0)       & (sdss['Z'] > sdss['Z_ERR'] ) # (sdss['SN_MEDIAN_ALL'] > 0.1 ) &

	error_reliable_boss = (boss[stellar_mass+'_up'] > boss[stellar_mass+'_low'] ) & (boss[stellar_mass+'_up'] > 0. ) & ( boss[stellar_mass+'_low'] > 0. ) & (boss[stellar_mass+'_up'] < 1e14 ) & ( boss[stellar_mass+'_low'] < 1e14 ) 
	error_reliable_sdss = (sdss[stellar_mass+'_up'] > sdss[stellar_mass+'_low'] ) & (sdss[stellar_mass+'_up'] > 0. ) & ( sdss[stellar_mass+'_low'] > 0. ) & (sdss[stellar_mass+'_up'] < 1e14 ) & ( sdss[stellar_mass+'_low'] < 1e14 ) 

	mass_reliable_boss = (boss[stellar_mass] > 1e6 ) & ( boss[stellar_mass] < 1e14 ) 
	mass_reliable_sdss = (sdss[stellar_mass] > 1e6 ) & ( sdss[stellar_mass] < 1e14 )
	
	ok_boss = (error_reliable_boss) & (mass_reliable_boss) & (redshift_reliable_boss)
	ok_sdss = (error_reliable_sdss) & (mass_reliable_sdss) & (redshift_reliable_sdss)

	snr = n.hstack(( boss['SN_MEDIAN_ALL'][ok_boss], sdss['SN_MEDIAN_ALL'][ok_sdss] ))
	zz = n.hstack(( boss['Z_NOQSO'][ok_boss], sdss['Z_NOQSO'][ok_sdss]))
	Ms = n.hstack(( n.log10(boss[stellar_mass][ok_boss]), n.log10(sdss[stellar_mass][ok_sdss]) ))
	err_moyenne = n.hstack(( abs((n.log10(boss[stellar_mass + '_up'][ok_boss]) - n.log10(boss[stellar_mass + '_low'][ok_boss]))/2.), abs((n.log10(sdss[stellar_mass + '_up'][ok_sdss]) - n.log10(sdss[stellar_mass + '_low'][ok_sdss]))/2.) ))
	out_s, outN_s = [], []
	for ii, sn in enumerate(sn_bins[:-1]):
		out = n.histogram(err_moyenne[(snr>sn_bins[ii])&(snr<sn_bins[ii+1])], bins = err_bins, normed=True)[0]
		outN = n.histogram(err_moyenne[(snr>sn_bins[ii])&(snr<sn_bins[ii+1])], bins = err_bins)[0]
		#N100 = (outN>1)
		out_s.append(out)
		outN_s.append(outN)

	return snr, zz, Ms, err_moyenne, sn_bins, x_err_bins, x_err_width, n.array(out_s), n.array(outN_s)


for imf in imfs:
	DATA = get_err_distribution(boss, sdss, imf = imf)
	f=open(os.path.join(out_dir, imf+"_pdf_DELTA_logM_SNR.pkl" ), 'w')
	cPickle.dump(DATA, f)
	f.close()
	# snr, zz, Ms, err_moyenne, sn_bins, x_err_bins, x_err_width, out_s, outN_s = DATA
	

for imf in imfs:
  f=open(os.path.join(out_dir, imf+"_pdf_DELTA_logM_SNR.pkl" ), 'r')
  snr, zz, Ms, err_moyenne, sn_bins, x_err_bins, x_err_width, out_s, outN_s = cPickle.load(f)
  f.close()
  p.figure(1, (4.5, 4.5))
  p.axes([0.2,0.2,0.7,0.7])
  cc=cycler('color', ['r', 'g', 'b', 'm', 'k'])
  for ii, (sn, c, out, outN) in enumerate(zip(sn_bins[:-1],cc, out_s, outN_s)):
    y=1.*n.array([n.sum(outN[:jj+1]) for jj in range(len(outN)-1)])/n.sum(outN)
    p.plot(x_err_bins[:-1], y, color=c['color'], ls='dashed')
    p.errorbar(x_err_bins[:-1], y, xerr=x_err_width[:-1], yerr=y*(y*outN[:1])**(-0.5), label=str(sn)+"-"+str(sn_bins[ii+1]), fmt='+', color=c['color'])
    
  p.legend(loc=0, frameon=False, title='SNR bins')
  p.ylabel('Cumulative normed distribution')
  p.xlim((-0.01, 1.5))
  p.ylim((-0.01, 1.01))
  p.grid()
  print("imf",imf)
  p.title(imf.split('_')[1])
  p.xlabel(r'$\sigma_{logM}$')
  p.savefig(os.path.join(out_dir, imf+"_pdf_DELTA_logM_SNR_cumulative.jpg" ))
  p.clf()


for degree in n.arange(1,5,1):
  for imf in imfs:
    f=open(os.path.join(out_dir, imf+"_pdf_DELTA_logM_SNR.pkl" ), 'r')
    snr, zz, Ms, err_moyenne, sn_bins, x_err_bins, x_err_width, out_s, outN_s = cPickle.load(f)
    f.close()

    p.figure(1, (4.5, 4.5))
    p.axes([0.2,0.2,0.7,0.7])
    cc=cycler('color', ['r', 'g', 'b', 'm', 'k'])
    print "DEGREE", degree
    for ii, (sn, c, out, outN) in enumerate(zip(sn_bins[:-1],cc, out_s, outN_s)):
      N100 = (outN>10)
      #print x_err_bins[N100], out[N100], x_err_width, out[N100]*outN[N100]**(-0.5)
      #p.plot(x_err_bins[N100], out[N100], marker='+', color=c['color'], label=str(sn_bins[ii])+"<SN<"+str(sn_bins[ii+1]))
      eb = p.errorbar(x_err_bins[N100], out[N100], xerr=x_err_width[N100], yerr = out[N100]*outN[N100]**(-0.5), label=str(sn)+"-"+str(sn_bins[ii+1]), fmt='+', color=c['color'])
      #p.fill_between(x_err_bins[N100], y1 = out[N100]-out[N100]*outN[N100]**(-0.5), y2 = out[N100]+out[N100]*outN[N100]**(-0.5), color=c['color'] ) 
      if len(out[N100])>3:
	print sn, "-", sn_bins[ii+1], " & " ,
	outP = n.polyfit(x_err_bins[N100], n.log10(out[N100]), deg=degree, w = 1/( outN[N100]**(-0.5)) )
	x_poly =n.arange(0.01,1.5,0.05)
	p.plot(x_poly, 10**n.polyval(outP, x_poly), c=c['color'], ls='dashed', lw=0.5)
	print n.round(outP,3) , " \\\\"
	#print x_err_bins[N100], out[N100]

    p.legend(loc=0, frameon=False, fontsize=12)
    p.yscale('log')
    p.ylabel('Normed distribution')
    p.xlim((-0.01, 1.5))
    p.ylim((1e-3, 10.))
    #p.title(imf+"degree="+str(degree))
    p.grid()
    p.xlabel(r'$\sigma_{logM}$')
    p.savefig(os.path.join(out_dir, imf+str(degree)+"_pdf_DELTA_logM_SNR.jpg" ))
    p.clf()


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




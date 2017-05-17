import os
import numpy as n
import astropy.io.fits as fits
import matplotlib
matplotlib.rcParams['agg.path.chunksize'] = 2000000
matplotlib.rcParams.update({'font.size': 13})
matplotlib.use('Agg')
import matplotlib.pyplot as p


"""
Plots log stellar mass vs. log(1+z) for each PROGRAMNAME.

"""


#prefix = 'SDSS'
#hdus = fits.open(os.path.join(os.environ['DATA_DIR'], 'spm', 'firefly', 'FireflyGalaxySdss26.fits'))
#redshift_reliable = (hdus[1].data['Z'] >= 0) & ( hdus[1].data['Z_ERR'] >= 0) & (hdus[1].data['ZWARNING'] == 0) & (hdus[1].data['Z'] > hdus[1].data['Z_ERR'] )

prefix = 'BOSS'
imf = 'Chabrier'
out_dir = os.path.join(os.environ['DATA_DIR'], 'spm', 'results', 'catalogs', imf)

hdus = fits.open(os.path.join(os.environ['DATA_DIR'], 'spm', 'firefly', 'FireflyGalaxyEbossDR14.fits'))
redshift_reliable =  (hdus[1].data['SN_MEDIAN_ALL'] > 0.1 ) & (hdus[1].data['CLASS_NOQSO'] == "GALAXY") & (hdus[1].data['Z_NOQSO'] >= 0) & ( hdus[1].data['Z_ERR_NOQSO'] >= 0) & (hdus[1].data['ZWARNING_NOQSO'] == 0) & (hdus[1].data['Z_NOQSO'] > hdus[1].data['Z_ERR_NOQSO'] )

def plot_all_prognames(hdus=hdus, imf=imf, prefix=prefix, out_dir = out_dir, redshift_reliable=redshift_reliable ) :
	stellar_mass = imf+'_stellar_mass'
	out_dir = os.path.join(os.environ['DATA_DIR'], 'spm', 'results', 'catalogs', imf)

	mass_reliable = (hdus[1].data[stellar_mass] > 0 ) & ( hdus[1].data[stellar_mass] < 13. ) & ( abs(hdus[1].data[stellar_mass + '_err_plus'] - hdus[1].data[stellar_mass]) < 0.4 ) & ( abs(hdus[1].data[stellar_mass + '_err_minus'] - hdus[1].data[stellar_mass]) < 0.4 )

	#good_plates = (hdus[1].data['PLATEQUALITY']=='good') &(hdus[1].data['TARGETTYPE']=='science')

	all_names = set(hdus[1].data['PROGRAMNAME'])
	all_names_arr = n.array(list(all_names))

	for ii in range(len(all_names_arr)):
		selection = (mass_reliable) & (redshift_reliable) & (hdus[1].data['PROGRAMNAME']==all_names_arr[ii])
		N_occ = len(selection.nonzero()[0])
		print all_names_arr[ii], N_occ
		if N_occ>1:
			p.figure(1, (4.5, 4.5))
			p.axes([0.2,0.2,0.7,0.7])
			p.plot(n.log10(1.+hdus[1].data['Z'][selection]), hdus[1].data[stellar_mass][selection], 'k+', rasterized=True, alpha=0.5) #, label=all_names_arr[ii]
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
			#p.plot((hdus[1].data['PLUG_RA'][selection]-180.)*n.pi/180., hdus[1].data['PLUG_DEC'][selection]*n.pi/180., 'k+', rasterized=True) # , label=all_names_arr[ii]
			p.plot(hdus[1].data['PLUG_RA'][selection], hdus[1].data['PLUG_DEC'][selection], 'k+', rasterized=True) # , label=all_names_arr[ii]
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
hdus = fits.open(os.path.join(os.environ['DATA_DIR'], 'spm', 'firefly', 'FireflyGalaxySdss26.fits'))
redshift_reliable = (hdus[1].data['SN_MEDIAN_ALL'] > 0.1 ) & (hdus[1].data['CLASS'] == "GALAXY") & (hdus[1].data['Z'] >= 0) & ( hdus[1].data['Z_ERR'] >= 0) & (hdus[1].data['ZWARNING'] == 0) & (hdus[1].data['Z'] > hdus[1].data['Z_ERR'] )

imf = 'Chabrier'
plot_all_prognames(hdus=hdus, imf=imf, prefix=prefix, redshift_reliable=redshift_reliable )
"""
imf = 'Salpeter'
plot_all_prognames(hdus=hdus, imf=imf, prefix=prefix, redshift_reliable=redshift_reliable )

imf = 'Kroupa'
plot_all_prognames(hdus=hdus, imf=imf, prefix=prefix, redshift_reliable=redshift_reliable )
"""




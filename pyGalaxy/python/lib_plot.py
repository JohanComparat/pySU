"""
This library enables to plot the different surveys on the same basis.
"""
import astropy.io.fits as fits
import numpy as n
import matplotlib
matplotlib.use('pdf')
matplotlib.rcParams['font.size']=14
import matplotlib.pyplot as p

def plotRaDecEBV(ra,dec,ebv,name,pDir):
	"""
	creates a figure with ra, dec, coded with E(B-V).
	:param ra: column right ascension
	:param dec: column declination
	:param ebv: column E(B-V)
	:param name: name of the figure.
	:param pDir: where the figure will be saved.
	"""
	p.figure(0,(6,5))
	p.axes([0.2,0.2,0.65,0.75])
	p.scatter(ra,dec,s=20,c=ebv,rasterized=True,edgecolor='none')
	p.xticks(rotation=70)
	p.xlabel('RA [deg, J2000]')
	p.ylabel('DEC [deg, J2000]')
	cb=p.colorbar()
	cb.set_label('E(B-V)')
	p.grid()
	p.savefig(pDir + name)
	p.clf()

def plotZ_SSR(zz,ssr,name,ylab,pDir):
	"""
	creates a figure with redshift, spectroscopic success rate.
	:param zz: column redshift
	:param ssr: column spectroscopic success rate
	:param name: name of the figure.
	:param ylab : y axis label
	:param pDir: where the figure will be saved.
	"""
	fig=p.figure(2,(5,5))
	p.axes([0.25,0.15,0.65,0.7])
	p.plot(zz,ssr,'b+',rasterized=True)
	p.xlabel('redshift')
	p.ylabel(ylab)
	#p.ylim((-17.5,-14.5))
	p.xlim((0.1,1.4))
	p.grid()
	p.savefig(pDir + name)
	p.clf()

def plot_I_TSR(mag,tsr,name,ylab,pDir):
	"""
	creates a figure with magnitude and target sampling rate.
	:param mag: column selection band magnitude
	:param tsr: column target sampling rate
	:param name: name of the figure.
	:param ylab : y axis label
	:param pDir: where the figure will be saved.
	"""
	fig=p.figure(2,(5,5))
	p.axes([0.25,0.15,0.65,0.7])
	p.plot(mag,tsr,'b+',rasterized=True)
	p.xlabel('selection magnitude')
	p.ylabel(ylab)
	#p.ylim((-17.5,-14.5))
	#p.xlim((0.1,1.4))
	p.grid()
	p.savefig(pDir + name)
	p.clf()



def plotZ_Flux(zz,ff,name,ylab,pDir):
	"""
	creates a figure with redshift, line flux.
	:param zz: column redshift
	:param ff: column line flux
	:param name: name of the figure.
	:param ylab : y axis label
	:param pDir: where the figure will be saved.
	"""
	bins=[n.arange(0.1,1.4,0.025),n.arange(-18,-14,0.1)]
	fig=p.figure(1,(6,5))
	p.axes([0.25,0.15,0.65,0.7])
	p.hist2d(zz,n.log10(ff),rasterized=True,bins=bins,cmin=2)
	p.xlabel('redshift')
	p.ylabel(ylab)
	cb=p.colorbar(shrink=0.8)
	cb.set_label('N ELG')
	p.ylim((-17.5,-14.5))
	p.xlim((0.1,1.4))
	p.grid()
	p.savefig(pDir + name)
	p.clf()

def plotZ_EW(zz,ew,name,ylab,pDir):
	"""
	creates a figure with redshift, equivalent width.
	:param zz: column redshift
	:param ew: column equivalent width
	:param name: name of the figure.
	:param ylab : y axis label
	:param pDir: where the figure will be saved.
	"""
	bins=[n.arange(0.1,1.4,0.025),n.arange(0,3,0.1)]
	fig=p.figure(1,(6,5))
	p.axes([0.25,0.15,0.65,0.7])
	p.hist2d(zz,n.log10(ew),rasterized=True,bins=bins,cmin=2)
	p.xlabel('redshift')
	p.ylabel(ylab)
	cb=p.colorbar(shrink=0.8)
	cb.set_label('N ELG')
	p.ylim((0,3))
	p.xlim((0.1,1.4))
	p.grid()
	p.savefig(pDir + name)
	p.clf()

def plotZ_Luminosity(zz,lum,name,ylab,pDir):
	"""
	creates a figure with redshift, luminosity.
	:param zz: column redshift
	:param ew: column luminosity
	:param name: name of the figure.
	:param ylab : y axis label
	:param pDir: where the figure will be saved.
	"""
	bins=[n.arange(0.1,1.4,0.025),n.arange(38,45,0.2)]
	fig=p.figure(1,(6,5))
	p.axes([0.25,0.15,0.65,0.7])
	p.hist2d(zz,n.log10(lum),rasterized=True,bins=bins,cmin=2)
	p.xlabel('redshift')
	p.ylabel(ylab)
	cb=p.colorbar(shrink=0.8)
	cb.set_label('N ELG')
	p.ylim((39,44))
	p.xlim((0.1,1.4))
	p.grid()
	p.savefig(pDir + name)
	p.clf()

def plot_EW_LF_measurement(lf_fits_file,lf_measurement_file,plotDir):
	"""
	plot the LF and histograms of EW and luminosity.
	"""
	zMean=lf_measurement_file.split('/')[-1][:-4].split('-')[-1][1:]
	line=lf_measurement_file.split('/')[-1][:-4].split('-')[0]
	lineDict = {'O2_3728' : r'$[O^{3728}_{II}]$', 'O3_5007' : r'$[O^{5007}_{III}]$', 'H1_4862' : r'$H^{4861}_{\beta}$', 'H1_6564': r'$H^{6564}_{\alpha}$'}
	lineLabel = lineDict[line]
	survey=lf_measurement_file.split('/')[-1][:-4].split('-')[1]
	hduG=fits.open(lf_fits_file)
	hdu=hduG[1]
	catalog=hdu.data
	completeness=hduG[0].header['COMPLETENESS']
	volume=hduG[0].header['VOLUME']
	#print completeness
	Lmin, Lmax, Lmean, phi, phi_err, phi_err_poisson, ngals= n.loadtxt( lf_measurement_file, unpack=True)
	# defines the main frame
	fig = p.figure(0,(7,7))
	fig.subplots_adjust(hspace=0.5,wspace=0.5)
	# defines the first panel
	fig.add_subplot(2,2,1)
	p.title(survey+" z="+zMean)
	tp=(phi>0)
	p.errorbar(Lmean[tp],phi[tp],yerr=phi_err[tp],xerr=[Lmean[tp]-Lmin[tp], Lmax[tp]-Lmean[tp]], fmt=None,elinewidth=1)
	p.axvline(completeness)
	p.axhline(10./volume)
	p.xlabel(r'$log_{10}(L$'+lineLabel+'$)$ [erg s$^{-1}$]')
	p.ylabel(r'$\Phi$ [Mpc$^{-3}/$dlog L]')
	p.yscale('log')
	p.xscale('log')
	p.xlim((1e39,1e44))
	p.ylim((1e-7,1e-1))
	p.grid()

	fig.add_subplot(2,2,2)
	bins= n.logspace(0,3.5,20)
	nn,bb,pp=p.hist(catalog[line+'_EW'],bins=bins,histtype='step')
	p.axvline(bb[n.argmax(nn)+1])
	p.yscale('log')
	p.xscale('log')
	p.xlim((1,1e4))
	p.ylim((0.8, len(catalog[line+'_EW'])/2. ))
	p.xlabel(r'$log_{10}(EW$'+lineLabel+'$)$ [A]')
	p.ylabel(r'counts')
	p.grid()

	fig.add_subplot(2,2,3)
	bins=n.logspace(39,44,20)
	p.hist(catalog[line+'_luminosity'],bins=bins,histtype='step')
	p.axvline(completeness)
	p.yscale('log')
	p.xscale('log')
	p.xlim((1e39,1e44))
	p.ylim((0.8, len(catalog[line+'_EW'])/2. ))
	p.xlabel(r'$log_{10}(L$'+lineLabel+'$)$ [A]')
	p.ylabel(r'counts')
	p.grid()

	fig.add_subplot(2,2,4)
	p.plot(catalog[line+'_EW'],catalog[line+'_luminosity'],'b+',rasterized=True)
	p.axhline(completeness)
	p.axvline(bb[n.argmax(nn)+1]*0.9,color='k',ls='dashed')
	p.axvline(bb[n.argmax(nn)+1],color='k')
	p.axvline(bb[n.argmax(nn)+1]*1.1,color='k',ls='dashed')
	p.xlabel(r'$log_{10}(EW$'+lineLabel+'$)$ [A]')
	p.ylabel(r'$log_{10}(L$'+lineLabel+'$)$ [A]')
	p.ylim((1e39,1e44))
	p.xlim((1,1e4))
	p.yscale('log')
	p.xscale('log')
	p.grid()

	p.savefig(plotDir+lf_measurement_file.split('/')[-1][:-4]+".pdf")
	p.clf()

def plot_EW_LF_measurement_simulation(lf_fits_file,lf_measurement_file,plotDir,lfsdata):
	"""
	plot the LF and histograms of EW and luminosity.
	"""
	zMean=lf_measurement_file.split('/')[-1][:-4].split('-')[-1][1:]
	line=lf_measurement_file.split('/')[-1][:-4].split('-')[0]
	lineDict = {'O2_3728' : r'$[O^{3728}_{II}]$', 'O3_5007' : r'$[O^{5007}_{III}]$', 'H1_4862' : r'$H^{4861}_{\beta}$', 'H1_6564': r'$H^{6564}_{\alpha}$'}
	lineLabel = lineDict[line]
	survey=lf_measurement_file.split('/')[-1][:-4].split('-')[1]
	hdu=fits.open(lf_fits_file)[1]
	catalog=hdu.data
	completeness=hdu.header['COMPLETENESS']
	#print completeness
	Lmin, Lmax, Lmean, phi, phi_err, ngals= n.loadtxt( lf_measurement_file, unpack=True)
	# defines the main frame
	fig = p.figure(0,(7,7))
	fig.subplots_adjust(hspace=0.5,wspace=0.5)
	# defines the first panel
	fig.add_subplot(2,2,1)
	p.title(survey+" z="+zMean)
	tp=(phi>0)
	p.errorbar(Lmean[tp],phi[tp],yerr=phi_err[tp], xerr=[Lmean[tp]-Lmin[tp], Lmax[tp]-Lmean[tp]], fmt=None,elinewidth=2,label="galform")
	#tp=(phi>0)
	#p.errorbar(Lmean[tp],phi[tp],yerr=phi_err_jackknife[tp],xerr=[Lmean[tp]-Lmin[tp], Lmax[tp]-Lmean[tp]], fmt=None,elinewidth=1)
	for lf in lfsdata:
		#print 'len lf', len(lf), lf
		if len(lf)==2:
			# case of a simulation points
			Lmin, Lmax, Lmean, phi, phi_err, ngals= n.loadtxt( lf[0], unpack=True)
			tp=(phi>0)
			p.errorbar(Lmean[tp],phi[tp],yerr=phi_err[tp], xerr=[Lmean[tp]-Lmin[tp], Lmax[tp]-Lmean[tp]], fmt=None,elinewidth=2,label=lf[1])

		if len(lf)==4:
			# case of data points
			Lmin, Lmax, Lmean, phi, phi_err, ngals= n.loadtxt( lf[0], unpack=True)
			p.errorbar(Lmean,phi*lf[2],yerr=phi*lf[2]*lf[3], xerr=[Lmean-Lmin, Lmax-Lmean], fmt=None,elinewidth=2,label=lf[1]+"Vcorr")
			Lmin, Lmax, Lmean, phi, phi_err, ngals= n.loadtxt( lf[0], unpack=True)
			tp=(phi>0)
			p.errorbar(Lmean[tp],phi[tp],yerr=phi_err_poisson[tp], xerr=[Lmean[tp]-Lmin[tp], Lmax[tp]-Lmean[tp]], fmt=None,elinewidth=2,label=lf[1])


	p.axvline(completeness)
	p.legend(fontsize=7,loc=3)
	p.xlabel(r'$log_{10}(L$'+lineLabel+'$)$ [erg s$^{-1}$]')
	p.ylabel(r'$\Phi$ [Mpc$^{-3}/$dlog L]')
	p.yscale('log')
	p.xscale('log')
	p.xlim((1e39,1e44))
	p.ylim((1e-7,1e-1))
	p.grid()

	fig.add_subplot(2,2,2)
	bins= n.logspace(0,3.5,20)
	nn,bb,pp=p.hist(catalog[line+'_EW'],bins=bins,histtype='step')
	p.axvline(bb[n.argmax(nn)+1])
	p.yscale('log')
	p.xscale('log')
	p.ylim((0.8, len(catalog[line+'_EW'])/2. ))
	p.xlabel(r'$log_{10}(EW$'+lineLabel+'$)$ [A]')
	p.ylabel(r'counts')
	p.grid()

	fig.add_subplot(2,2,3)
	bins=n.logspace(39,44,20)
	p.hist(catalog[line+'_luminosity'],bins=bins,histtype='step')
	p.axvline(completeness)
	p.yscale('log')
	p.xscale('log')
	p.xlim((1e39,1e44))
	p.ylim((0.8, len(catalog[line+'_EW'])/2. ))
	p.xlabel(r'$log_{10}(L$'+lineLabel+'$)$ [A]')
	p.ylabel(r'counts')
	p.grid()

	fig.add_subplot(2,2,4)
	p.plot(catalog[line+'_EW'],catalog[line+'_luminosity'],'b+',rasterized=True)
	p.axhline(completeness)
	p.axvline(bb[n.argmax(nn)])
	p.xlabel(r'$log_{10}(EW$'+lineLabel+'$)$ [A]')
	p.ylabel(r'$log_{10}(L$'+lineLabel+'$)$ [A]')
	p.ylim((1e39,1e44))
	p.xlim((1e-3,1e4))
	p.yscale('log')
	p.xscale('log')
	p.grid()

	p.savefig(plotDir+lf_measurement_file.split('/')[-1][:-4]+".pdf")
	p.clf()

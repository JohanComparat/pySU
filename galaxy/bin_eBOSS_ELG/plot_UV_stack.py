#! /usr/bin/env python

"""
This script produces the stacks for emission line luminosity limited samples.
"""
import sys
import os 
from os.path import join
import glob
import numpy as n
import astropy.io.fits as fits
import SpectraStackingEBOSS as sse
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

line_list_abs = n.array([2249.88, 2260.78, 2344.21, 2374.46, 2382.76, 2576.88, 2586.65, 2594.50, 2600.17, 2606.46, 2796.35, 2803.53, 2852.96])
line_list_abs_names = n.array(['FeII', 'FeII', 'FeII', 'FeII', 'FeII', 'MnII', 'FeII', 'MnII','FeII', 'MnII', 'MgII','MgII','MgI'])
line_list_em = n.array([2327, 2365.55, 2396.36, 2612.65,2626.45])
line_list_em_names = n.array(['CII]', 'FeII*', 'FeII*', 'FeII*', 'FeII*'])

def plot_me(qty):
	# create all stacks
	dataList_UV = n.array(glob.glob(join(os.environ['HOME'], "SDSS/stacks", "eboss-elg_*_"+qty+"_*.UVstack")))
	dataList = n.array(glob.glob(join(os.environ['HOME'], "SDSS/stacks", "eboss-elg_*_"+qty+"_*.stack")))
	
	fig=p.figure(0,(7.2, 13.7), frameon=False)
	
	fig.add_subplot(311, xlim=((2240, 2410)))
	for specList in dataList_UV:
		bn = os.path.basename(specList)[10:-8].split('_')
		print(bn)
		if qty == 'fast_lmass':
			bnl = str(n.round(float(bn[0]),3))+'<z<'+str(n.round(float(bn[2]),3))+', '+str(n.round(float(bn[3]),3))+'<'+qty+'<'+str(n.round(float(bn[6]),3))
		else:
			bnl = str(n.round(float(bn[0]),3))+'<z<'+str(n.round(float(bn[2]),3))+', '+str(n.round(float(bn[3]),3))+'<'+qty+'<'+str(n.round(float(bn[5]),3))
		dd=fits.open(specList)[1].data
		wl=dd['wavelength'         ]
		s1 = (wl>2200)&(wl<2900)&(dd['NspectraPerPixel'   ]>0.5*n.max(dd['NspectraPerPixel'   ]))
		p.plot(dd['wavelength'         ][s1], dd['medianStack'        ][s1], label= bnl, lw=0.7 )

	for xx, nn in zip(line_list_abs, line_list_abs_names ):
		p.plot([xx,xx],[0,1], ls='dashed', color='k')
		p.text(xx,0.7,nn,rotation=90)

	for xx, nn in zip(line_list_em, line_list_em_names ):
		p.axvline([xx,xx],[1,2], ls='dashed', color='g')
		p.text(xx,1.1,nn,rotation=90, color='g')

	p.legend(frameon=False)
	p.grid()

	fig.add_subplot(312, ylabel=r'Flux/Fcont]', xlim=((2570, 2640)))
	for specList in dataList_UV:
		dd=fits.open(specList)[1].data
		wl=dd['wavelength'         ]
		s1 = (wl>2200)&(wl<2900)&(dd['NspectraPerPixel'   ]>0.5*n.max(dd['NspectraPerPixel'   ]))
		p.plot(dd['wavelength'         ][s1], dd['medianStack'        ][s1], lw=0.7 )

	for xx, nn in zip(line_list_abs, line_list_abs_names ):
		p.axvline([xx,xx],[0,1], ls='dashed', color='k')
		p.text(xx,0.7,nn,rotation=90)

	for xx, nn in zip(line_list_em, line_list_em_names ):
		p.axvline([xx,xx],[1,2], ls='dashed', color='g')
		p.text(xx,1.1,nn,rotation=90, color='g')
	p.grid()


	fig.add_subplot(313, xlabel='wavelength [Angstrom, rest frame]', xlim=((2780, 2870)))
	for specList in dataList_UV:
		dd=fits.open(specList)[1].data
		wl=dd['wavelength'         ]
		s1 = (wl>2200)&(wl<2900)&(dd['NspectraPerPixel'   ]>0.5*n.max(dd['NspectraPerPixel'   ]))
		p.plot(dd['wavelength'         ][s1], dd['medianStack'        ][s1], lw=0.7 )

	for xx, nn in zip(line_list_abs, line_list_abs_names ):
		p.axvline([xx,xx],[0,1], ls='dashed', color='k')
		p.text(xx,0.7,nn,rotation=90)

	for xx, nn in zip(line_list_em, line_list_em_names ):
		p.axvline([xx,xx],[1,2], ls='dashed', color='g')
		p.text(xx,1.1,nn,rotation=90, color='g')
	
	p.grid()

	p.tight_layout()
	p.savefig(join(os.environ['HOME'], "SDSS/stacks", "eboss-elg_"+qty+".UVstack")+".png")
	p.clf()


	p.figure(0,(9,5))
	p.title(qty)
	for specList in dataList:
		bn = os.path.basename(specList)[10:-8]
		dd=fits.open(specList)[1].data
		wl=dd['wavelength'         ]
		s1 = (dd['NspectraPerPixel'   ]>0.5*n.max(dd['NspectraPerPixel'   ]))
		p.plot(dd['wavelength'         ][s1], dd['medianStack'        ][s1], label= bn )
		#dd['meanStack'          ][s1]
		#dd['meanWeightedStack'  ][s1]
		#dd['jackknifeSpectra'   ][s1]
		#dd['jackknifStackErrors'][s1]
		#dd['NspectraPerPixel'   ][s1]

	p.legend(frameon=False)
	p.tight_layout()
	p.savefig(join(os.environ['HOME'], "SDSS/stacks", "eboss-elg_"+qty+".stack")+".png")
	p.clf()

plot_me(qty = 'fast_lmass' )
plot_me(qty = 'g'          )
plot_me(qty = 'gr'         )
plot_me(qty = 'rz'         )

os.system("cp -r ~/SDSS/stacks/*.png ~/wwwDir/sdss/elg/stacks/")

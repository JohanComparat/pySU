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
	
	p.figure(0,(9,5))
	p.title(qty)
	for specList in dataList_UV:
		bn = os.path.basename(specList)[10:-8]
		dd=fits.open(specList)[1].data
		wl=dd['wavelength'         ]
		s1 = (wl>2200)&(wl<2900)&(dd['NspectraPerPixel'   ]>0.5*n.max(dd['NspectraPerPixel'   ]))
		p.plot(dd['wavelength'         ][s1], dd['medianStack'        ][s1], label= bn )
		#dd['meanStack'          ][s1]
		#dd['meanWeightedStack'  ][s1]
		#dd['jackknifeSpectra'   ][s1]
		#dd['jackknifStackErrors'][s1]
		#dd['NspectraPerPixel'   ][s1]


	for xx, nn in zip(line_list_abs, line_list_abs_names ):
		p.axvline(xx, ls='dashed', color='k')
		p.text(xx,0.3,nn,rotation=90)

	for xx, nn in zip(line_list_em, line_list_em_names ):
		p.axvline(xx, ls='dashed', color='=g')
		p.text(xx,1.1,nn,rotation=90, color='g')


	p.legend(frameon=False)
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

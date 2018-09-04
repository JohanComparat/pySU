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

qty = 'fast_lmass'
qty = 'g'
qty = 'gr'
qty = 'rz'


# create all stacks
dataList_UV = n.array(glob.glob(join(os.environ['HOME'], "SDSS/lss/catalogs/3", "eboss-elg_*_"+qty+"_*.UVstack")))
dataList = n.array(glob.glob(join(os.environ['HOME'], "SDSS/lss/catalogs/3", "eboss-elg_*_"+qty+"_*.stack")))
 
p.figure(0,(9,5))
p.title(qty)
for specList in dataList_UV:
	bn = os.path.basename(specList)[10:-8]
	dd=fits.open(specList)[1].data
	wl=dd['wavelength'         ]
	s1 = (wl>2000)(wl<3700)&(	dd['NspectraPerPixel'   ]>0.5*n.max(dd['NspectraPerPixel'   ]))
	p.plot(dd['wavelength'         ][s1], dd['medianStack'        ][s1], label= bn )
	#dd['meanStack'          ][s1]
	#dd['meanWeightedStack'  ][s1]
	#dd['jackknifeSpectra'   ][s1]
	#dd['jackknifStackErrors'][s1]
	#dd['NspectraPerPixel'   ][s1]

p.legend(frameon=False)
p.tight_layout()
p.savefig(specList".png")
p.clf()


p.figure(0,(9,5))
p.title(qty)
for specList in dataList_UV:
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
p.savefig(specList".png")
p.clf()
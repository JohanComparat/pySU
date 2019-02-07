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
from scipy.integrate import quad

# ICM GCM lines absorption

line_names_ICM_CGM_abs = {
	'Mg_1_2853':  'MgI 2853',
	'Mg_2_2796': 'MgII 2796',
	'Mg_2_2804': 'MgII 2804',
	'Mn_2_2577': 'MgII 2577', 
	'Mn_2_2594': 'MgII 2594', 
	'Mn_2_2606': 'MgII 2606',
    'Fe_2_2600': 'FeII 2600',
	'Fe_2_2587': 'FeII 2587',
	'Fe_2_2383': 'FeII 2383',
	'Fe_2_2374': 'FeII 2374',
	'Fe_2_2344': 'FeII 2344',
	'Fe_2_2261': 'FeII 2261',
	'Fe_2_2250': 'FeII 2250'
}

line_lambda_ICM_CGM_abs = {
	'Mg_1_2853': 2852.96, 
	'Mg_2_2796': 2796.35,
	'Mg_2_2804': 2803.53,
	'Mn_2_2577': 2576.88, 
	'Mn_2_2594': 2594.50, 
	'Mn_2_2606': 2606.46,
	'Fe_2_2600': 2600.17, 
	'Fe_2_2587': 2586.55, 
	'Fe_2_2383': 2382.76, 
	'Fe_2_2374': 2374.46, 
	'Fe_2_2344': 2344.21, 
	'Fe_2_2261': 2260.78,
	'Fe_2_2250': 2249.88
}

# ICM GCM lines emission
line_names_ICM_CGM_abs = {
	'Fes_2_2632': 'FeII* 2632',
	'Fes_2_2626': 'FeII* 2626',
	'Fes_2_2613': 'FeII* 2613',
	'Fes_2_2396': 'FeII* 2396',
	'Fes_2_2381': 'FeII* 2381',
	'Fes_2_2366': 'FeII* 2366',
	'Fes_2_2281': 'FeII* 2281',
	'Fes_2_2270': 'FeII* 2270',
	
}

line_lambda_ICM_CGM_abs = {
	'Fes_2_2632': 2632.11,
	'Fes_2_2626': 2626.45,
	'Fes_2_2613': 2612.65,
	'Fes_2_2396': 2396.35,
	'Fes_2_2381': 2381.49,
	'Fes_2_2366': 2365.55,
	'Fes_2_2281': 2280.62,
	'Fes_2_2270': 2269.52,
}

# nebular emission
line_names_nebular_emi = {
	'C2_2324': 'CII] 2324',
	'C2_2325': 'CII] 2325',
	'C2_2326': 'CII] 2326',
	'C2_2328': 'CII] 2328',
	'C2_2329': 'CII] 2329',
	'Ne4_2423': '[NeIV] 2423', 
	'Ne4_2425': '[NeIV] 2425', 
	'O2_2470': '[OII] 2470', 
	'O2_2471': '[OII] 2471' 
}

line_lambda_nebular_emi = {
	'C2_2324': 2324.21, 
	'C2_2325': 2325.40, 
	'C2_2326': 2326.11, 
	'C2_2328': 2327.64, 
	'C2_2329': 2328.83, 
	'Ne4_2423': 2422.56, 
	'Ne4_2425': 2425.14,
	'O2_2470': 2470.97, 
	'O2_2471': 2471.09 
	}

def normalize_C0(x,y,yerr):
	"""
	We then mask out
	absorption and emission features and fit a cubic polyno-
	mial function through the rest of the spectrum. 
	Using
	the best-fit polynomial function as an estimate of the un-
	derlying continuum, F lambda
	we normalize the observed spectrum to obtain the continuum-normalized spectrum
	"""
	# masking sky contaminated pixels
	maskLambda = n.loadtxt(os.path.join(os.environ['GIT_SPM'],'data',"dr12-sky-mask.txt"), unpack=True)
	ratio = n.min(abs(10000.*n.log10(n.outer(x, 1./maskLambda))), axis=1)
	margin = 1.5
	veto_sky = ( ratio <= margin )
	
	# UV mask
	UV_mask = (x>2000)&(x<3600)
	
	# UV line mask
	ratio = n.min(abs(10000.*n.log10(n.outer(x, 1./line_list_abs))), axis=1)
	margin = 10
	veto_line_abs = ( ratio <= margin )

	ratio = n.min(abs(10000.*n.log10(n.outer(x, 1./line_list_em))), axis=1)
	margin = 10
	veto_line_em = ( ratio <= margin )
	
	# MASKING BAD DATA
	bad_data = n.isnan(y) | n.isinf(y) | (y <= 0.0) | n.isnan(yerr) | n.isinf(yerr)
	# creating new arrays
	x = x[(UV_mask)&(veto_sky==False)&(bad_data==False)&(veto_line_abs==False)&(veto_line_em==False)] 
	y = y[(UV_mask)&(veto_sky==False)&(bad_data==False)&(veto_line_abs==False)&(veto_line_em==False)] 
	yerr = yerr[(UV_mask)&(veto_sky==False)&(bad_data==False)&(veto_line_abs==False)&(veto_line_em==False)] 
	
	out=n.polyfit(x, y, 4, w=1/yerr)
	return out

def W0_abs(x, y, l_min, l_max):
	it = interp1d(x,1-y)
	return quad(it, l_min, l_max)[0]

def W0_emi(x, y, l_min, l_max):
	it = interp1d(x,y-1)
	return quad(it, l_min, l_max)[0]


def open_spec(path_2_file):
	return fits.open(path_2_file)[1].data

# list of spectra
dataList_UV = n.array(glob.glob(join(os.environ['HOME'], "SDSS/stacks", "eboss-elg_*_"+qty+"_*.UVstack")))

def create_xy_data(path_2_file):
	# opens spectrum
	dd=open_spec(path_2_file)
	# normalizes with a fourth order polynomial
	wl=dd['wavelength'         ]
	s1 = (wl>2200)&(wl<2900)&(dd['NspectraPerPixel'   ]>0.5*n.max(dd['NspectraPerPixel'   ]))
	#, label= bnl, lw=0.7 )
	pfit = normalize_C0(dd['wavelength'         ][s1], dd['medianStack'        ][s1] , dd['jackknifStackErrors'][s1])
	Fcont = n.polyval(pfit, dd['wavelength'         ][s1])
	wave = dd['wavelength'         ][s1]
	spec = dd['medianStack'        ][s1] / Fcont
	specErr = dd['jackknifStackErrors'        ][s1] / Fcont
	return wave, spec, specErr

idx=0
x, y, ye = create_xy_data(dataList_UV[idx])
# now fitting line features :


def plot_v_profile(x,y,ye):
	# Figure 8.
	fig=p.figure(0, (5., 8.), frameon=False)
	# top panel
	fig.add_subplot(111, ylabel=r'F/Fcont', xlim=((-800, 600)), ylim=((0.9,2)))
	#p.plot( x-line_lambda_ICM_CGM_abs['Fes_2_2632'],y,label=line_names_ICM_CGM_abs['Fes_2_2632'] )
	p.plot( x-line_lambda_ICM_CGM_abs['Fes_2_2626'],y,label=line_names_ICM_CGM_abs['Fes_2_2626'] )
	p.plot( x-line_lambda_ICM_CGM_abs['Fes_2_2613'],y,label=line_names_ICM_CGM_abs['Fes_2_2613'] )
	p.plot( x-line_lambda_ICM_CGM_abs['Fes_2_2396'],y,label=line_names_ICM_CGM_abs['Fes_2_2396'] )
	p.plot( x-line_lambda_ICM_CGM_abs['Fes_2_2381'],y,label=line_names_ICM_CGM_abs['Fes_2_2381'] )
	p.plot( x-line_lambda_ICM_CGM_abs['Fes_2_2366'],y,label=line_names_ICM_CGM_abs['Fes_2_2366'] )
	p.plot( x-line_lambda_ICM_CGM_abs['Fes_2_2281'],y,label=line_names_ICM_CGM_abs['Fes_2_2281'] )
	p.plot( x-line_lambda_ICM_CGM_abs['Fes_2_2270'],y,label=line_names_ICM_CGM_abs['Fes_2_2270'] )
	p.legend(frameon=False)
	p.grid()

	fig.add_subplot(112, ylabel=r'F/Fcont', xlim=((-800, 600)), ylim=((0.9,2)))
	
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

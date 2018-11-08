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


agn_file = join(os.environ['HOME'],"SDSS/stacks/v1", "xagn_0.38.stack")
agn=fits.open(agn_file)[1].data
s1 = (agn['wavelength'         ]>2200)&(agn['wavelength'         ]<2900)
xA, yA = agn['wavelength'         ][s1], agn['medianStack'        ][s1]

agn=fits.open("/home/comparat/data-4most/opsim_data/survey_team_inputs/Xgal/templates/4most_AGN_type2_zmin_00_zmax_03_EBV_0_01.fits")[1].data

s1 = (agn['LAMBDA']/(1.15)>2200)&(agn['LAMBDA']/(1.15)<3800)
xA, yA = agn['LAMBDA'][s1]/(1.15), agn['FLUX'][s1]/(3.5950675388804168e-16)

ELG=fits.open("/home/comparat/stacks/zhu_2015/eBOSS_ELG_NUV_composite.fits")[1].data
ELG_a=fits.open("/home/comparat/stacks/zhu_2015/eBOSS_ELG_composite.fits")[1].data
#eBOSS_ELG_composite.fits  


line_list_abs = n.array([2249.88, 2260.78, 2344.21, 2374.46, 2382.76, 2576.88, 2586.65, 2594.50, 2600.17, 2606.46, 2796.35, 2803.53, 2852.96])
line_list_abs_names = n.array(['FeII', 'FeII', 'FeII', 'FeII', 'FeII', 'MnII', 'FeII', 'MnII','FeII', 'MnII', 'MgII','MgII','MgI'])
line_list_em = n.array([2327, 2365.55, 2396.36, 2612.65,2626.45])
line_list_em_names = n.array(['CII]', 'FeII*', 'FeII*', 'FeII*', 'FeII*'])

def plot_me(qty='O2EW'):
	# create all stacks
	# O2EW
	# O2lum
	dataList_UV = n.array(glob.glob(join(os.environ['HOME'], "stacks", "eboss-elg_0.6_z_1.2_*_"+qty+"_*.UVstack")))
	dataList = n.array([el[:-7]+"stack" for el in dataList_UV])
	print(dataList_UV,dataList)
	fig=p.figure(0,(7.2, 13.7), frameon=False)
	
	fig.add_subplot(311, xlim=((2240, 2410)), ylim=((0,2)))
	p.plot(xA,yA,'k',lw=0.5, label='NLAGN')
	for specList in dataList_UV:
		bn = os.path.basename(specList)[10:-8].split('_')
		print(bn)
		#if qty == 'fast_lmass':
			#bnl = str(n.round(float(bn[0]),3))+'<z<'+str(n.round(float(bn[2]),3))+', '+str(n.round(float(bn[3]),3))+'<'+qty+'<'+str(n.round(float(bn[6]),3))
		#else:
		bnl = str(n.round(float(bn[0]),3))+'<z<'+str(n.round(float(bn[2]),3))+', '+str(n.round(float(bn[3]),3))+'<'+qty+'<'+str(n.round(float(bn[5]),3))
		dd=fits.open(specList)[1].data
		wl=dd['wavelength']
		s1 = (wl>2200)&(wl<2900)&(dd['NspectraPerPixel'   ]>0.5*n.max(dd['NspectraPerPixel'   ]))
		#if float(bn[0])>0.5:
		p.plot(dd['wavelength'         ][s1], dd['medianStack'        ][s1], label= bnl, lw=0.7 )

	for xx, nn in zip(line_list_abs, line_list_abs_names ):
		p.plot(n.array([xx,xx]),n.array([0,1]), ls='dashed', color='k', lw=0.5)
		p.text(xx,0.2,nn,rotation=90)

	for xx, nn in zip(line_list_em, line_list_em_names ):
		p.plot(n.array([xx,xx]),n.array([1,2]), ls='dashed', color='g', lw=0.5)
		p.text(xx,1.4,nn,rotation=90, color='g')
	
	ok = (ELG['WAVE']>2200)&(ELG['WAVE']<2900)
	p.plot(ELG['WAVE'][ok], ELG['FLUXMEDIAN'][ok] , label='Zhu15')

	p.legend(frameon=False)
	p.grid()

	fig.add_subplot(312, ylabel=r'F/Fcont', xlim=((2570, 2640)), ylim=((0,2)))
	p.plot(xA,yA,'k',lw=0.5, label='NLAGN')
	for specList in dataList_UV:
		dd=fits.open(specList)[1].data
		wl=dd['wavelength'         ]
		s1 = (wl>2200)&(wl<2900)&(dd['NspectraPerPixel'   ]>0.5*n.max(dd['NspectraPerPixel'   ]))
		bn = os.path.basename(specList)[10:-8].split('_')
		#if float(bn[0])>0.5:
		p.plot(dd['wavelength'         ][s1], dd['medianStack'        ][s1], lw=0.7 )

	for xx, nn in zip(line_list_abs, line_list_abs_names ):
		p.plot(n.array([xx,xx]),n.array([0,1]), ls='dashed', color='k', lw=0.5)
		p.text(xx,0.2,nn,rotation=90)

	for xx, nn in zip(line_list_em, line_list_em_names ):
		p.plot(n.array([xx,xx]),n.array([1,2]), ls='dashed', color='g', lw=0.5)
		p.text(xx,1.4,nn,rotation=90, color='g')
	p.grid()
	ok = (ELG['WAVE']>2200)&(ELG['WAVE']<2900)
	p.plot(ELG['WAVE'][ok], ELG['FLUXMEDIAN'][ok] , label='Zhu15')


	fig.add_subplot(313, xlabel='wavelength [Angstrom, rest frame]', xlim=((2780, 2870)), ylim=((0,2)))
	p.plot(xA,yA,'k',lw=0.5, label='NLAGN')
	for specList in dataList_UV:
		dd=fits.open(specList)[1].data
		wl=dd['wavelength'         ]
		s1 = (wl>2200)&(wl<2900)&(dd['NspectraPerPixel'   ]>0.5*n.max(dd['NspectraPerPixel'   ]))
		bn = os.path.basename(specList)[10:-8].split('_')
		#if float(bn[0])>0.5:
		p.plot(dd['wavelength'         ][s1], dd['medianStack'        ][s1], lw=0.7 )

	for xx, nn in zip(line_list_abs, line_list_abs_names ):
		p.plot(n.array([xx,xx]),n.array([0,1]), ls='dashed', color='k', lw=0.5)
		p.text(xx,0.2,nn,rotation=90)

	for xx, nn in zip(line_list_em, line_list_em_names ):
		p.plot(n.array([xx,xx]),n.array([1,2]), ls='dashed', color='g', lw=0.5)
		p.text(xx,1.4,nn,rotation=90, color='g')
	
	p.grid()
	ok = (ELG['WAVE']>2200)&(ELG['WAVE']<2900)
	p.plot(ELG['WAVE'][ok], ELG['FLUXMEDIAN'][ok] , label='Zhu15')

	p.tight_layout()
	p.savefig(join(os.environ['HOME'], "stacks", "eboss-elg_"+qty+".UVstack")+".png")
	p.clf()

	p.figure(1,(9,5))
	p.title(qty)
	for specList in dataList:
		bn = os.path.basename(specList)[10:-8].split('_')
		print(bn)
		#if qty == 'fast_lmass':
			#bnl = str(n.round(float(bn[0]),3))+'<z<'+str(n.round(float(bn[2]),3))+', '+str(n.round(float(bn[3]),3))+'<'+qty+'<'+str(n.round(float(bn[6]),3))
		#else:
		bnl = str(n.round(float(bn[0]),3))+'<z<'+str(n.round(float(bn[2]),3))+', '+str(n.round(float(bn[3]),3))+'<'+qty+'<'+str(n.round(float(bn[5]),3))
		dd=fits.open(specList)[1].data
		wl=dd['wavelength'         ]
		s1 = (dd['NspectraPerPixel'   ]>0.5*n.max(dd['NspectraPerPixel'   ]))
		p.plot(dd['wavelength'         ][s1], dd['medianStack'        ][s1], label= bnl, lw=0.5 )
		#dd['meanStack'          ][s1]
		#dd['meanWeightedStack'  ][s1]
		#dd['jackknifeSpectra'   ][s1]
		#dd['jackknifStackErrors'][s1]
		#dd['NspectraPerPixel'   ][s1]

	p.plot(ELG['WAVE'], ELG['FLUXMEDIAN'], label='Zhu15')
	p.grid()
	p.xlim((2300,3800))
	p.plot(xA,yA,'k',lw=0.5, label='NLAGN')
	p.yscale('log')
	p.legend(frameon=False)
	p.tight_layout()
	p.savefig(join(os.environ['HOME'], "stacks", "eboss-elg_"+qty+".stack")+".png")
	p.clf()
	
plot_me(qty = 'O2EW' )
plot_me(qty = 'O2lum' )

#plot_me(qty = 'mass' )
#plot_me(qty = 'g'          )
#plot_me(qty = 'gr'         )
#plot_me(qty = 'rz'         )
#plot_me(qty = 'rw1'         )
os.system("rm ~/wwwDir/sdss/elg/stacks/*.png")
os.system("cp -r ~/stacks/*.png ~/wwwDir/sdss/elg/stacks/")

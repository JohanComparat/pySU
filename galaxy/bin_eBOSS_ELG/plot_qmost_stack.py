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
import lineListVac as ll

em_line_list = [
	[1240.14, 'N V' ,  'darkgreen'],
	[1305.53, 'O I' ,  'darkgreen'],
	[1335.31, 'C II', 'darkgreen' ],
	[1397.61, 'Si IV', 'darkgreen' ],
	[1399.8, 'Si IV + O IV', 'darkgreen' ],
	[ll.C4_1548, r'C IV', 'darkgreen'],
	[1640.42, 'He II', 'darkgreen'],
	[1750.26, 'N III]', 'darkgreen'],
	[ll.C3_1908  , r'C III', 'darkgreen' ],
	[2327.00, 'CII]', 'darkgreen'],
	[2396.36, 'FeII*', 'darkgreen'],
	[2626.45, 'FeII*', 'darkgreen'],
	[3346.82, '[Ne V]', 'darkgreen'],
	[3426.84, '[Ne V]', 'darkgreen'],
	[ll.O2_mean  , r'[O II]', 'darkgreen'],
	[3759.99, '[Fe VII]', 'darkgreen'],
	[ll.Ne3_3869 , r'[Ne III]', 'darkgreen'],
    # [ll.Ne3_3968 , r'[Ne III]', 'darkgreen'],
    [ll.O3_4363  , r'[O III]'  , 'darkgreen'],
    [ll.O3_4960  , r'[O III]'  , 'darkgreen'],
    [ll.O3_5007  , r'[O III]'  , 'darkgreen'],
	[5160.33, '[Fe VII]', 'darkgreen'],
	[ll.O1_5578  , r'O I', 'darkgreen' ],
	[5722.30, '[Fe VII]', 'darkgreen'],
	[5877.29, 'He I', 'darkgreen'],
	[6087.98, '[Fe VII]', 'darkgreen'],
	[ll.O1_6302  , r'O I' , 'darkgreen'],
    [ll.O1_6365  , r'O I' , 'darkgreen'],
    [ll.N2_5756  , r'[N II]' , 'darkgreen'],
    [ll.N2_6549  , r'[N II]' , 'darkgreen'],
    [ll.N2_6585  , r'[N II]' , 'darkgreen'],
    [ll.S2_6718  , r'[S II]', 'darkgreen'],
    [ll.S2_6732  , r'[S II]', 'darkgreen'],
    [ll.Ar3_7137 , r'[Ar III]' , 'darkgreen'],
    ]

abs_line_list = [
	[911.753, r'Ly$_{limit}$', 'black'],
	[1025.7220, r'Ly$_\beta$', 'black'],
	[ll.H1_1216, r'Ly$_\alpha$', 'black'],
	[1857.40, 'Al III', 'darkgreen'],
	#
	[2344.21, 'FeII', 'darkgreen'],
	[2382.76, 'Fe II', 'darkgreen'],
	[2600.17, 'FeII', 'darkgreen'],
	[2798.75, 'MgII', 'darkgreen'],
	#
	[3835.397, r'H$\eta$', 'black'],
	[3889.064, r'H$\zeta$', 'black'],
	[3934.777, 'Ca(K)', 'magenta'],
	[3969.588, 'Ca(H)', 'magenta'],
	[ll.H1_3970 , r'H$_\epsilon$', 'black'],
	#
    [ll.H1_4102 , r'H$_\delta$', 'black'],
	[4305.61, 'G', 'magenta'],
    [ll.H1_4341 , r'H$_\gamma$', 'black'],
    [ll.He2_4686 , r'He II', 'darkgreen'],
    [ll.H1_4862 , r'H$_\beta$', 'black'],
	#
	[5176.7, 'MgI b', 'magenta'],
	[ll.He2_5411, r'He II', 'darkgreen'],
	[5895.6, r'NaI D$_{1,2}$', 'magenta'],
	[ll.H1_6564 , r'H$_\alpha$', 'black'],
	#
	[8500.36, 'Ca II', 'magenta'],
	[8544.44, 'Ca II', 'magenta'],
	[8664.52, 'Ca II', 'magenta'],
]

# line_list_abs = n.array([      2249.88, 2260.78, 2344.21, 2374.46, 2382.76, 2576.88, 2586.65, 2594.50, 2600.17, 2606.46, 2796.35, 2803.53, 2852.96])
# line_list_abs_names = n.array(['FeII' ,  'FeII',  'FeII',  'FeII',  'FeII',  'MnII',  'FeII',  'MnII',  'FeII',  'MnII',  'MgII',  'MgII',   'MgI'])
# line_list_em = n.array([2327, 2365.55, 2396.36, 2612.65,2626.45])
# line_list_em_names = n.array(['CII]', 'FeII*', 'FeII*', 'FeII*', 'FeII*'])

#stack_dir = join( os.environ['HOME'], "SDSS/stacks/v2" )
stack_dir = join( os.environ['HOME'], "SDSS/stacks" )
file_out =  join(stack_dir,"X_AGN", "DR16_ELG-stitched-stack.fits")
file_out =  join(stack_dir,"X_AGN", "ROSAT_AGNT1-stitched-stack.fits")
file_out =  join(stack_dir,"X_AGN", "ROSAT_AGNT2-stitched-stack.fits")

def plot_spec( p_2_stack = file_out ):
	print('plots', p_2_stack)
	fig=p.figure(7, (14.0, 14.0), frameon=False)
	fig.add_subplot(411, ylabel=r'F$_\lambda$')
	stack = fits.open(p_2_stack)[1].data
	s1 = (stack['wavelength']>0)
	stack = stack[s1]
	y_min = n.min(stack['medianStack'])
	y_max = n.max(stack['medianStack'])
	delta_y = y_max - y_min
	p.xlim((n.min(stack['wavelength']), n.max(stack['wavelength'])))
	p.ylim((y_min - delta_y * 0.2 , y_max + delta_y * 0.2 ))
	# p.yscale('log')
	# lines above
	for elem in em_line_list:
		print(elem)
		if elem[0]>n.min(stack['wavelength'][5]) and elem[0]<n.max(stack['wavelength'][-5]) :
			xpos = n.searchsorted(stack['wavelength'], elem[0])
			ypos = n.max(stack['medianStack'][xpos-10:xpos+10]) + delta_y * 0.1
			# p.plot(n.array([elem[0], elem[0]]), em_dash_Y, ls='dashed', color='k', lw=0.5)
			p.text(elem[0], ypos, r'$^{----}$' + elem[1], rotation=90, c='darkgreen')
	# lines below
	for elem in abs_line_list:
		print(elem)
		if elem[0]>n.min(stack['wavelength'][5]) and elem[0]<n.max(stack['wavelength'][-5]) :
			xpos = n.searchsorted(stack['wavelength'], elem[0])
			ypos = n.min(stack['medianStack'][xpos-30:xpos+30]) - delta_y * 0.2
			# p.plot(n.array([elem[0], elem[0]]), em_dash_Y, ls='dashed', color='k', lw=0.5)
			p.text(elem[0], ypos, elem[1] + r'$^{---}$', rotation=90, c='magenta')
	p.plot(stack['wavelength'], stack['medianStack'], lw=0.7)
	p.grid()
	p.tight_layout()
	#
	print('standard deviation')
	fig.add_subplot(412, ylabel=r'per cent')
	stack = fits.open(p_2_stack)[1].data
	s1 = (stack['wavelength']>0)
	stack = stack[s1]
	y_min = n.min( [ stack['jackknifStackErrors'], stack['NspectraPerPixel']**-0.5 ] )
	y_max = n.max( [ stack['jackknifStackErrors'], stack['NspectraPerPixel']**-0.5 ] )
	p.xlim((n.min(stack['wavelength']), n.max(stack['wavelength'])))
	p.ylim(( y_min/1.1 , y_max*1.1 ))
	p.plot(stack['wavelength'], stack['jackknifStackErrors']/stack['medianStack'], lw=0.7, label=r'$\sigma^{var}_{JK}$')
	p.plot(stack['wavelength'], stack['NspectraPerPixel']**-0.5, lw=2, label=r'$1/\sqrt{N}$')
	p.grid()
	p.legend()
	p.yscale('log')
	p.tight_layout()
	print('correlation coefficient')
	fig.add_subplot(212, ylabel='Wavelength rest-frame [Angstrom]',  xlabel='Wavelength rest-frame [Angstrom]')
	CR = n.corrcoef(stack['jackknifeSpectra'])
	WLa = n.array([ stack['wavelength'] for el in stack['wavelength'] ])
	WLb = WLa.T
	highCorr_sel = ( abs(CR) > 0.8 ) & (CR>0)
	xx = WLa[highCorr_sel]
	yy = WLb[highCorr_sel]
	cr_val = CR[highCorr_sel]
	p.scatter(xx, yy, c=cr_val, s=1, rasterized = True)
	p.colorbar(shrink=0.8)
	p.tight_layout()
	p.savefig(p_2_stack+".png")
	p.clf()

plot_spec(  join(stack_dir,"X_AGN", "ROSAT_AGNT1-stitched-stack.fits") )

plot_spec(  join(stack_dir,"X_AGN", "ROSAT_AGNT2-stitched-stack.fits") )
plot_spec(  join(stack_dir,"X_AGN", "DR16_ELG-stitched-stack.fits") )


sys.exit()

p.figure(1,(9,5))
p.title(qty)
#p.plot(xA,yA/n.median(yA),'k',lw=0.5, label='Narrow Line AGN')
ok = (ELG_a['WAVE']>2300)&(ELG_a['WAVE']<3800)

p.plot(ELG_a['WAVE'][ok], ELG_a['FLUXMEDIAN'][ok]/ELG_a['FLUXMEDIAN_ERR'][ok], lw=0.3, label='Zhu15 SNR pilot study')
for specList in dataList_UV:
	bn = os.path.basename(specList)[10:-8].split('_')
	bnl = str(n.round(float(bn[0]),3))+'<z<'+str(n.round(float(bn[2]),3))+', '+str(n.round(float(bn[3]),3))+'<'+qty+'<'+str(n.round(float(bn[5]),3))
	dd=fits.open(specList)[1].data
	wl=dd['wavelength']
	s1 = (dd['NspectraPerPixel'   ]>0.5*n.max(dd['NspectraPerPixel'   ]))
	p.plot(dd['wavelength'][s1], 200*dd['medianStack'        ][s1]/n.median(dd['medianStack'        ][s1]), label= bnl )
	#dd['meanStack'          ][s1]
	#dd['meanWeightedStack'  ][s1]
	#dd['jackknifeSpectra'   ][s1]
	#dd['jackknifStackErrors'][s1]
	#dd['NspectraPerPixel'   ][s1]

p.plot(dd['wavelength'][s1], dd['medianStack'][s1]/dd['jackknifStackErrors'][s1], lw=0.3, label='SNR full eBOSS')
p.grid()
p.xlim((2300,3800))
p.yscale('log')
p.ylabel('S/N per pixel')
p.xlabel('wavelength [A]')
p.legend(frameon=False)
p.tight_layout()
p.savefig(p_2_stack + ".png")
p.clf()


p.figure(1,(9,5))
p.title(qty)
#p.plot(xA,yA/n.median(yA),'k',lw=0.5, label='Narrow Line AGN')
for specList in dataList_UV:
	bn = os.path.basename(specList)[10:-8].split('_')
	bnl = str(n.round(float(bn[0]),3))+'<z<'+str(n.round(float(bn[2]),3))+', '+str(n.round(float(bn[3]),3))+'<'+qty+'<'+str(n.round(float(bn[5]),3))
	dd=fits.open(specList)[1].data
	wl=dd['wavelength']
	s1 = (dd['NspectraPerPixel'   ]>0.5*n.max(dd['NspectraPerPixel'   ]))
	p.plot(dd['wavelength'][s1], dd['medianStack_UVnormed'        ][s1], label= bnl )
	#dd['meanStack'          ][s1]
	#dd['meanWeightedStack'  ][s1]
	#dd['jackknifeSpectra'   ][s1]
	#dd['jackknifStackErrors'][s1]
	#dd['NspectraPerPixel'   ][s1]

#ok = (ELG_a['WAVE']>2300)&(ELG_a['WAVE']<3800)
#p.plot(ELG_a['WAVE'][ok], ELG_a['FLUXMEDIAN'][ok]/n.median(ELG_a['FLUXMEDIAN'][ok]), label='Zhu15')
p.grid()
p.xlim((2240,2400))
p.ylim((0,2))
p.legend(frameon=False)
p.tight_layout()
p.savefig(join(stack_dir, "eboss-elg-2240-2400_"+qty+".stack")+".png")
p.clf()


p.figure(1,(9,5))
p.title(qty)
#p.plot(xA,yA/n.median(yA),'k',lw=0.5, label='Narrow Line AGN')
for specList in dataList_UV:
	bn = os.path.basename(specList)[10:-8].split('_')
	bnl = str(n.round(float(bn[0]),3))+'<z<'+str(n.round(float(bn[2]),3))+', '+str(n.round(float(bn[3]),3))+'<'+qty+'<'+str(n.round(float(bn[5]),3))
	dd=fits.open(specList)[1].data
	wl=dd['wavelength']
	s1 = (dd['NspectraPerPixel'   ]>0.5*n.max(dd['NspectraPerPixel'   ]))
	p.plot(dd['wavelength'][s1], dd['medianStack_UVnormed'        ][s1], label= bnl )
	#dd['meanStack'          ][s1]
	#dd['meanWeightedStack'  ][s1]
	#dd['jackknifeSpectra'   ][s1]
	#dd['jackknifStackErrors'][s1]
	#dd['NspectraPerPixel'   ][s1]

#ok = (ELG_a['WAVE']>2300)&(ELG_a['WAVE']<3800)
#p.plot(ELG_a['WAVE'][ok], ELG_a['FLUXMEDIAN'][ok]/n.median(ELG_a['FLUXMEDIAN'][ok]), label='Zhu15')
p.grid()
p.xlim((2570,2640))
p.ylim((0,2))
p.legend(frameon=False)
p.tight_layout()
p.savefig(join(stack_dir, "eboss-elg-2570-2640_"+qty+".stack")+".png")
p.clf()


p.figure(1,(9,5))
p.title(qty)
#p.plot(xA,yA/n.median(yA),'k',lw=0.5, label='Narrow Line AGN')
for specList in dataList_UV:
	bn = os.path.basename(specList)[10:-8].split('_')
	bnl = str(n.round(float(bn[0]),3))+'<z<'+str(n.round(float(bn[2]),3))+', '+str(n.round(float(bn[3]),3))+'<'+qty+'<'+str(n.round(float(bn[5]),3))
	dd=fits.open(specList)[1].data
	wl=dd['wavelength']
	s1 = (dd['NspectraPerPixel'   ]>0.5*n.max(dd['NspectraPerPixel'   ]))
	p.plot(dd['wavelength'][s1], dd['medianStack_UVnormed'        ][s1], label= bnl )
	#dd['meanStack'          ][s1]
	#dd['meanWeightedStack'  ][s1]
	#dd['jackknifeSpectra'   ][s1]
	#dd['jackknifStackErrors'][s1]
	#dd['NspectraPerPixel'   ][s1]

#ok = (ELG_a['WAVE']>2300)&(ELG_a['WAVE']<3800)
#p.plot(ELG_a['WAVE'][ok], ELG_a['FLUXMEDIAN'][ok]/n.median(ELG_a['FLUXMEDIAN'][ok]), label='Zhu15')
p.grid()
p.xlim((2780,2870))
p.ylim((0,2))
p.legend(frameon=False)
p.tight_layout()
p.savefig(join(stack_dir, "eboss-elg-2780-2870_"+qty+".stack")+".png")
p.clf()

plot_me(qty = 'O2EW' )
#plot_me(qty = 'O2lum' )

#plot_me(qty = 'mass' )
#plot_me(qty = 'g'          )
#plot_me(qty = 'gr'         )
#plot_me(qty = 'rz'         )
#plot_me(qty = 'rw1'         )
os.system("rm ~/wwwDir/sdss/elg/stacks/*.png")
os.system("cp -r "+stack_dir+"/*.png ~/wwwDir/sdss/elg/stacks/")

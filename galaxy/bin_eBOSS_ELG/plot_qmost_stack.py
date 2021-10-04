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
	[1240.14, 'N V'],
	[ll.C4_1548, r'C IV'],
	[1640.42, 'He II'],
	[1750.26, 'N III]'],
	[ll.C3_1908  , r'C III'   ],
	[2327.00, 'CII]'],
	[2396.36, 'FeII*'],
	[2626.45, 'FeII*'],
	[3346.82, '[Ne V]'],
	[3426.84, '[Ne V]'],
	[ll.O2_mean  , r'[O II]'  ],
	[3759.99, '[Fe VII]'],
	[ll.Ne3_3869 , r'[Ne III]'],
    # [ll.Ne3_3968 , r'[Ne III]'],
    [ll.O3_4363  , r'[O III]'  ],
    [ll.O3_4960  , r'[O III]'  ],
    [ll.O3_5007  , r'[O III]'  ],
	[5160.33, '[Fe VII]'],
	[ll.O1_5578  , r'O I' ],
	[5722.30, '[Fe VII]'],
	[5877.29, 'He I'],
	[6087.98, '[Fe VII]'],
	[ll.O1_6302  , r'O I' ],
    [ll.O1_6365  , r'O I' ],
    [ll.N2_5756  , r'[N II]' ],
    [ll.N2_6549  , r'[N II]' ],
    [ll.N2_6585  , r'[N II]' ],
    [ll.S2_6718  , r'[S II]'],
    [ll.S2_6732  , r'[S II]'],
    [ll.Ar3_7137 , r'[Ar III]' ],
    ]

abs_line_list = [
	[1857.40, 'Al III'],
	[2344.21, 'FeII'],
	[2382.76, 'Fe II'],
	[2600.17, 'FeII'],
	[2798.75, 'MgII'],
	[3933.70, 'Ca(H)'],
	[3968.50, 'Ca(K)'],
	[ll.H1_3970 , r'H$_\epsilon$'],
    [ll.H1_4102 , r'H$_\delta$'],
	[4304.4, 'G band'],
    [ll.H1_4341 , r'H$_\gamma$'],
    [ll.H1_4862 , r'H$_\beta$'],
	[5175.30, 'Mg'],
	[5894.00, 'Na'],
	[ll.H1_6564 , r'H$_\alpha$'],
    [ll.H1_1216 , r'Ly$_\alpha$'],
    [ll.He2_4686 , r'He II'],
    [ll.He2_5411 , r'He II'],
    ]

line_list_abs = n.array([      2249.88, 2260.78, 2344.21, 2374.46, 2382.76, 2576.88, 2586.65, 2594.50, 2600.17, 2606.46, 2796.35, 2803.53, 2852.96])
line_list_abs_names = n.array(['FeII' ,  'FeII',  'FeII',  'FeII',  'FeII',  'MnII',  'FeII',  'MnII',  'FeII',  'MnII',  'MgII',  'MgII',   'MgI'])
line_list_em = n.array([2327, 2365.55, 2396.36, 2612.65,2626.45])
line_list_em_names = n.array(['CII]', 'FeII*', 'FeII*', 'FeII*', 'FeII*'])

#stack_dir = join( os.environ['HOME'], "SDSS/stacks/v2" )
stack_dir = join( os.environ['HOME'], "SDSS/stacks" )
file_out =  join(stack_dir,"X_AGN", "DR16_ELG-stitched-stack.fits")
file_out =  join(stack_dir,"X_AGN", "ROSAT_AGNT1-stitched-stack.fits")
file_out =  join(stack_dir,"X_AGN", "ROSAT_AGNT2-stitched-stack.fits")

def plot_me( p_2_stack = file_out ):
	print('plots', p_2_stack)
	fig=p.figure(2, (14.0, 7.0), frameon=False)
	fig.add_subplot(111, ylabel=r'F$_\lambda$')#, xlim=((2240, 2410)), ylim=((0.7,1.2)))

	stack = fits.open(p_2_stack)[1].data
	s1 = (stack['wavelength']>0)
	stack = stack[s1]

	y_min = n.min(stack['medianStack'])
	y_max = n.max(stack['medianStack'])
	delta_y = y_max - y_min

	# t_height_em = y_max + delta_y * 0.1 # * 1.2
	# em_dash_Y = n.array([y_max+ delta_y * 0.05, t_height_em])
	#
	# t_height_abs = y_min * 0.75
	# ab_dash_Y = n.array([t_height_abs, y_min,])

	p.xlim((n.min(stack['wavelength']), n.max(stack['wavelength'])))
	p.ylim((y_min - delta_y * 0.2 , y_max + delta_y * 0.2 ))
	# p.yscale('log')

	# usual line list
	for elem in em_line_list:
		print(elem)
		if elem[0]>n.min(stack['wavelength'][5]) and elem[0]<n.max(stack['wavelength'][-5]) :
			xpos = n.searchsorted(stack['wavelength'], elem[0])
			ypos = n.max(stack['medianStack'][xpos-10:xpos+10]) + delta_y * 0.1
			# p.plot(n.array([elem[0], elem[0]]), em_dash_Y, ls='dashed', color='k', lw=0.5)
			p.text(elem[0], ypos, r'$^{----}$' + elem[1], rotation=90, c='darkgreen')

	# usual line list
	for elem in abs_line_list:
		print(elem)
		if elem[0]>n.min(stack['wavelength'][5]) and elem[0]<n.max(stack['wavelength'][-5]) :
			xpos = n.searchsorted(stack['wavelength'], elem[0])
			ypos = n.min(stack['medianStack'][xpos-30:xpos+30]) - delta_y * 0.2
			# p.plot(n.array([elem[0], elem[0]]), em_dash_Y, ls='dashed', color='k', lw=0.5)
			p.text(elem[0], ypos, elem[1] + r'$^{---}$', rotation=90, c='magenta')

	# # emission list in the UV
	# for xx, nn in zip(line_list_em, line_list_em_names ):
	# 	p.plot(n.array([xx,xx]), em_dash_Y, ls='dashed', color='g')
	# 	p.text(xx,t_height_em,nn,rotation=90, color='g')
	# # absorption list in the UV
	# for xx, nn in zip(line_list_abs, line_list_abs_names ):
	# 	p.plot(n.array([xx,xx]), ab_dash_Y, ls='dashed', color='k')
	# 	p.text(xx, t_height_abs, nn,rotation=90)
	p.plot(stack['wavelength'], stack['medianStack'], lw=0.7)

	# p.legend(frameon=False)
	p.grid()
	#
	# xmin, xmax = 2100, 3000
	# fig.add_subplot(312, ylabel=r'F$_\lambda$', xlim=((xmin, xmax )) )
	#
	# stack = fits.open(p_2_stack)[1].data
	# s1 = ( stack['wavelength'] > xmin ) &  ( stack['wavelength'] < xmax )
	# stack = stack[s1]
	#
	# y_min = n.min(stack['medianStack'])
	# y_max = n.max(stack['medianStack'])
	#
	# t_height_em = y_max * 1.2
	# em_dash_Y = n.array([y_max*0.8, t_height_em])
	#
	# t_height_abs = y_min * 0.75
	# ab_dash_Y = n.array([t_height_abs, y_min,])
	#
	# p.xlim((n.min(stack['wavelength']), n.max(stack['wavelength'])))
	# p.ylim((y_min * 0.7, y_max * 1.5))
	# # p.yscale('log')
	#
	# # usual line list
	# for elem in line_list:
	# 	p.plot(n.array([elem[0], elem[0]]), em_dash_Y, ls='dashed', color='k')
	# 	p.text(elem[0], t_height_em, elem[1], rotation=90)
	#
	# # emission list in the UV
	# for xx, nn in zip(line_list_em, line_list_em_names ):
	# 	p.plot(n.array([xx,xx]), em_dash_Y, ls='dashed', color='g')
	# 	p.text(xx,t_height_em,nn,rotation=90, color='g')
	# # absorption list in the UV
	# for xx, nn in zip(line_list_abs, line_list_abs_names ):
	# 	p.plot(n.array([xx,xx]), ab_dash_Y, ls='dashed', color='k')
	# 	p.text(xx, t_height_abs, nn,rotation=90)
	# p.plot(stack['wavelength'], stack['medianStack'], lw=0.7)
	#
	# # p.legend(frameon=False)
	# p.grid()
	#
	#
	# xmin, xmax = 3700, 5100
	# fig.add_subplot(313, ylabel=r'F$_\lambda$', xlabel='wavelength [Angstrom, rest frame]', xlim=((xmin, xmax )) )
	#
	# stack = fits.open(p_2_stack)[1].data
	# s1 = ( stack['wavelength'] > xmin ) &  ( stack['wavelength'] < xmax )
	# stack = stack[s1]
	#
	# y_min = n.min(stack['medianStack'])
	# y_max = n.max(stack['medianStack'])
	#
	# t_height_em = y_max * 1.2
	# em_dash_Y = n.array([y_max*0.8, t_height_em])
	#
	# t_height_abs = y_min * 0.75
	# ab_dash_Y = n.array([t_height_abs, y_min,])
	#
	# p.xlim((n.min(stack['wavelength']), n.max(stack['wavelength'])))
	# p.ylim((y_min * 0.7, y_max * 1.5))
	# # p.yscale('log')
	#
	# # usual line list
	# for elem in line_list:
	# 	p.plot(n.array([elem[0], elem[0]]), em_dash_Y, ls='dashed', color='k')
	# 	p.text(elem[0], t_height_em, elem[1], rotation=90)
	#
	# p.plot(stack['wavelength'], stack['medianStack'], lw=0.7)
	#
	# # p.legend(frameon=False)
	# p.grid()
	p.tight_layout()
	p.savefig(p_2_stack+".png")
	p.clf()

plot_me(  join(stack_dir,"X_AGN", "DR16_ELG-stitched-stack.fits") )
plot_me(  join(stack_dir,"X_AGN", "ROSAT_AGNT1-stitched-stack.fits") )
plot_me(  join(stack_dir,"X_AGN", "ROSAT_AGNT2-stitched-stack.fits") )


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

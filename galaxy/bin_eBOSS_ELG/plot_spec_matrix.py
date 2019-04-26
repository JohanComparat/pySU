import numpy as n
import glob
import os, sys
import astropy.io.fits as fits
import time
t0=time.time()
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p
import numpy as n

from cycler import cycler
# 1. Setting prop cycle on default rc parameter
p.rc('lines', linewidth=1.3)
p.rc('axes', prop_cycle=cycler('color', ["#638bd9", "#e586b6", "#dd7c25"]) )#, "#598664", "#9631a9", "#cff159", "#534f55", "#ab1519", "#89dbde"]) )

path_2_fly_dir = os.path.join(os.environ['HOME'], 'wwwDir/sdss/elg/')


spec_dir = os.path.join(os.environ['HOME'],"SDSS/stacks")
#specList = os.path.join(spec_dir, "eboss-elg_0.2_z_1.5.asc") 
out_file = os.path.join(spec_dir, 'eboss-elg_0.2_z_1.5'+".specMatrix")
for IDX_j in n.arange(0, 253454, 4096):
	print(IDX_j, time.time()-t0)

def resample_it(IDX_j):
	#IDX_j = 122880
	IDX_min=IDX_j
	IDX_max=IDX_j+4096
	IDX_str=str(IDX_min).zfill(6)+'-'+str(IDX_max).zfill(6)
	FLUXES     = n.loadtxt(out_file+'.'+IDX_str+'.dat')
	med_val = n.median(FLUXES, axis=1)
	F_med = FLUXES/med_val
	out = n.median(F_med.reshape(128,32,4096), axis=1)
	wavelength = n.loadtxt(out_file+'.wavelength.'+IDX_str+'.dat')
	data       = n.loadtxt(out_file+'.shapes.'+IDX_str+'.dat')
	plates, mjds, fiberids, redshifts = n.loadtxt(out_file+'.list.'+IDX_str+'.dat')
	out_z = n.median(redshifts.reshape(128, 32), axis=1)
	return out, out_z, wavelength


path_2_figure = os.path.join( path_2_fly_dir, 'matrixSpecPlot.resamp32.png')

out_0, out_z_0, wavelength = resample_it(IDX_j)
out_1, out_z_1, wavelength = resample_it(IDX_j+4096)
out_2, out_z_2, wavelength = resample_it(IDX_j+4096*2)
out_3, out_z_3, wavelength = resample_it(IDX_j+4096*3)

out = n.vstack(( out_0, out_1, out_2, out_3 ))
out_z = n.hstack(( out_z_0, out_z_1, out_z_2, out_z_3 ))

#xx = n.searchsorted(wavelength, n.array([2800., 3728., 5007., ])*(1+redshifts[0]))
#yy = n.zeros(len(xx))
#tt = n.array(['MgII 2800', '[OII] 3728', '[OIII] 5007' ])
#out[out>100.]=100.
#out[out<0.01]=0.01
fig = p.figure(2, (14.5, 6.5) )
#p.imshow(n.log10(FLUXES/med_val), cmap='tab20', rasterized = True)
p.imshow(n.log10(out), cmap='gist_gray', rasterized = True, vmin=-1, vmax=2.)
p.colorbar(shrink=0.7)
idx = n.arange(0, 4096, 512)
p.xticks(idx, wavelength[idx].astype('int'))
idy = n.arange(0, 512, 64)
p.yticks(idy, n.round(out_z[idy],2))
p.xlabel('wavelength')
p.ylabel('redshift')
#for x,y,t in zip(xx,yy,tt):
	#p.plot(x,y,'ro', markersize=5)
	#p.text(x,y-10,t,color='b', fontsize=14, rotation=0)
p.savefig(path_2_figure)
p.clf()

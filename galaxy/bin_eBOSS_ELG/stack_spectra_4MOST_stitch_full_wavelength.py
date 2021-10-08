#! /usr/bin/env python

"""
This script produces the stacks for samples defined by a list of spectra identifiers.

nohup python3 stack_spectra_ELG.py > stack_spectra_ELG.log &

cd ~/software/linux/qmost_templates/data
cp ~/SDSS/stacks/*/*stitched-stack.fits .
# push to the git
"""
import sys
import os 
from os.path import join
import glob
import numpy as n
import SpectraStackingEBOSS as sse
import astropy.io.fits as fits
from astropy.table import Table, Column

# # T2 templates
# ~/SDSS/stacks/stack_dir$ ls -lh ROSAT_AGNT2*stack
# ROSAT_AGNT2_zmin_00_zmax_02..stack for 6100 AA < lambda < 8100 AA
# ROSAT_AGNT2_zmin_00_zmax_05..stack for 3650 AA < lambda < 6100 AA
# ROSAT_AGNT2_zmin_03_zmax_08..stack for 3000 AA < lambda < 3650 AA
# ROSAT_AGNT2_zmin_05_zmax_10..stack for 2300 AA < lambda < 3000 AA
#
# adjust heights of each template so they match in the boundary region
#
# Write a template in the 4MOST format.
#
# # T1 template
# ~/SDSS/stacks/X_AGN
# ROSAT_AGNT1_zmin_*..stack

stack_dir = join(os.environ['HOME'],"SDSS/stacks")



###########
###########
# LRG DR16
###########
###########

T_stacks = n.array( glob.glob( join(stack_dir, "*LRG*stack") ) )
T_stacks.sort()
cww_e = '/home/comparat/data1/templates/lephare_dev/sed/GAL/CWW_KINNEY/CWW_E_ext.sed'
T_stacks = n.hstack(( T_stacks, '/home/comparat/SDSS/stacks/X_AGN/ROSAT_clusterGAL_zmin_03_zmax_08..stack', cww_e ))
file_out =  join(stack_dir,"X_AGN", "DR16LRG-stitched-stack.fits")

frac_Nspec = 0.001

idx_list = [9, 0, 7]
wl_mins = [7600, 3400, 2500 ]
wl_maxs = [12000, 7600, 3400 ]

ii=-1
jj = idx_list[ii]
wl0, wl1  = wl_mins[ii], wl_maxs[ii]
print(T_stacks[jj], wl0, wl1)
hdu = fits.open(T_stacks[jj])[1].data
s_wl = (hdu['wavelength'] > wl0) & (hdu['wavelength'] < wl1) & (hdu['NspectraPerPixel'] > n.max(hdu['NspectraPerPixel']) * frac_Nspec)
med_norm = n.median(hdu['medianStack'][s_wl][-10:])
t0 = Table(hdu[s_wl])
t0['medianStack'] = t0['medianStack']/med_norm
t0['meanStack'] = t0['meanStack']/med_norm
t0['jackknifStackErrors'] = t0['jackknifStackErrors']/med_norm
t0['jackknifeSpectra'] = t0['jackknifeSpectra']/med_norm
t0['NspectraPerPixelMax'] = n.max(t0['NspectraPerPixel'])* n.ones_like(t0['NspectraPerPixel'])

ts = [ t0 ]
for jj, wl0, wl1 in zip(idx_list[::-1][1:], wl_mins[::-1][1:], wl_maxs[::-1][1:]):
	print(T_stacks[jj], wl0, wl1)
	if T_stacks[jj][-3:]=='sed':
		data_i = n.loadtxt(T_stacks[jj], unpack = True)
		s_wl = (data_i[0]>wl0 ) & ( data_i[1] < wl1)
		data = data_i.T[s_wl].T
		med_norm = n.median(data[1][:10] )
		previous_norm = n.median( ts[-1]['medianStack'][-10:] )
		m_factor = previous_norm/med_norm
		t1 = Table()
		t1['wavelength'] = data[0]
		t1['medianStack'] = data[1]*m_factor
		t1['meanStack'] = data[1]*m_factor
		t1['jackknifStackErrors'] = -99. * n.ones( len(data[1]) )
		t1['jackknifeSpectra']    = -99. * n.ones( ( len(data[1]), 10 ) )
		t1['NspectraPerPixel']    = -99. * n.ones( len(data[1]) )
		t1['NspectraPerPixelMax'] = -99. * n.ones( len(data[1]) )
	else:
		hdu = fits.open(T_stacks[jj])[1].data
		s_wl = (hdu['wavelength']>wl0 ) & ( hdu['wavelength'] < wl1) & (hdu['NspectraPerPixel']>n.max(hdu['NspectraPerPixel']) * frac_Nspec )
		med_norm = n.median( hdu['medianStack'][s_wl][:10] )
		previous_norm = n.median( ts[-1]['medianStack'][-10:] )
		m_factor = previous_norm/med_norm
		t1 = Table(hdu[s_wl])
		t1['medianStack'] = t1['medianStack'] * m_factor
		t1['meanStack'] = t1['meanStack'] * m_factor
		t1['jackknifStackErrors'] = t1['jackknifStackErrors'] * m_factor
		t1['jackknifeSpectra'] = t1['jackknifeSpectra'] * m_factor
		t1['NspectraPerPixelMax'] = n.max(t1['NspectraPerPixel']) * n.ones_like(t1['NspectraPerPixel'])
	ts.append(t1)

t_out = Table( n.hstack(ts) )
t_out.remove_column('meanStack')
t_out.write(file_out, overwrite = True)

###########
###########
# QSO DR16
###########
###########

T_stacks = n.array( glob.glob( join(stack_dir,"X_AGN", "ROSAT_AGNT1*stack") ) )
T_stacks.sort()
T_stacks = n.hstack(( T_stacks, '/home/comparat/SDSS/stacks/STACK_DR16_QSO_zmin_35_zmax_59.ascii.stack' ))
file_out =  join(stack_dir,"X_AGN", "ROSAT_AGNT1-DR16QSO-stitched-stack.fits")

frac_Nspec = 0.5

idx_list = [0, 1, 11, 21, 31, 39, 40]
wl_mins = [6900, 3200, 2500, 1800, 1480, 1350, 950 ]
wl_maxs = [8100, 6900, 3200, 2500, 1800, 1480, 1350]

ii=-1
jj = idx_list[ii]
wl0, wl1  = wl_mins[ii], wl_maxs[ii]
print(T_stacks[jj], wl0, wl1)
hdu = fits.open(T_stacks[jj])[1].data
s_wl = (hdu['wavelength'] > wl0) & (hdu['wavelength'] < wl1) & (hdu['NspectraPerPixel'] > n.max(hdu['NspectraPerPixel']) * frac_Nspec)
med_norm = n.median(hdu['medianStack'][s_wl][-10:])
t0 = Table(hdu[s_wl])
t0['medianStack'] = t0['medianStack']/med_norm
t0['meanStack'] = t0['meanStack']/med_norm
t0['jackknifStackErrors'] = t0['jackknifStackErrors']/med_norm
t0['jackknifeSpectra'] = t0['jackknifeSpectra']/med_norm
t0['NspectraPerPixelMax'] = n.max(t0['NspectraPerPixel'])* n.ones_like(t0['NspectraPerPixel'])

ts = [ t0 ]
for jj, wl0, wl1 in zip(idx_list[::-1][1:], wl_mins[::-1][1:], wl_maxs[::-1][1:]):
	print(T_stacks[jj], wl0, wl1)
	hdu = fits.open(T_stacks[jj])[1].data
	s_wl = (hdu['wavelength']>wl0 ) & ( hdu['wavelength'] < wl1) & (hdu['NspectraPerPixel']>n.max(hdu['NspectraPerPixel']) * frac_Nspec )
	med_norm = n.median( hdu['medianStack'][s_wl][:10] )
	previous_norm = n.median( ts[-1]['medianStack'][-10:] )
	m_factor = previous_norm/med_norm
	t1 = Table(hdu[s_wl])
	t1['medianStack'] = t1['medianStack'] * m_factor
	t1['meanStack'] = t1['meanStack'] * m_factor
	t1['jackknifStackErrors'] = t1['jackknifStackErrors'] * m_factor
	t1['jackknifeSpectra'] = t1['jackknifeSpectra'] * m_factor
	t1['NspectraPerPixelMax'] = n.max(t1['NspectraPerPixel']) * n.ones_like(t1['NspectraPerPixel'])
	ts.append(t1)

t_out = Table( n.hstack(ts) )
t_out.remove_column('meanStack')
t_out.write(file_out, overwrite = True)

###########
###########
# ELG DR16
###########
###########

T_stacks = n.array( glob.glob( join(stack_dir, "STACK_DR16_ELG*stack") ) )
T_stacks.sort()
file_out =  join(stack_dir,"X_AGN", "DR16_ELG-stitched-stack.fits")

frac_Nspec = 0.5

idx_list = [0, 2, 7, 1]
wl_mins = [5300, 2950, 2200, 1750 ]
wl_maxs = [8100, 5300, 2950, 2200 ]

ii=-1
jj = idx_list[ii]
wl0, wl1  = wl_mins[ii], wl_maxs[ii]
print(T_stacks[jj], wl0, wl1)
hdu = fits.open(T_stacks[jj])[1].data
s_wl = (hdu['wavelength'] > wl0) & (hdu['wavelength'] < wl1) & (hdu['NspectraPerPixel'] > n.max(hdu['NspectraPerPixel']) * frac_Nspec)
med_norm = n.median(hdu['medianStack'][s_wl][-10:])
t0 = Table(hdu[s_wl])
t0['medianStack'] = t0['medianStack']/med_norm
t0['meanStack'] = t0['meanStack']/med_norm
t0['jackknifStackErrors'] = t0['jackknifStackErrors']/med_norm
t0['jackknifeSpectra'] = t0['jackknifeSpectra']/med_norm
t0['NspectraPerPixelMax'] = n.max(t0['NspectraPerPixel'])* n.ones_like(t0['NspectraPerPixel'])

ts = [ t0 ]
for jj, wl0, wl1 in zip(idx_list[::-1][1:], wl_mins[::-1][1:], wl_maxs[::-1][1:]):
	print(T_stacks[jj], wl0, wl1)
	hdu = fits.open(T_stacks[jj])[1].data
	s_wl = (hdu['wavelength']>wl0 ) & ( hdu['wavelength'] < wl1) & (hdu['NspectraPerPixel']>n.max(hdu['NspectraPerPixel']) * frac_Nspec )
	med_norm = n.median( hdu['medianStack'][s_wl][:10] )
	previous_norm = n.median( ts[-1]['medianStack'][-10:] )
	m_factor = previous_norm/med_norm
	t1 = Table(hdu[s_wl])
	t1['medianStack'] = t1['medianStack'] * m_factor
	t1['meanStack'] = t1['meanStack'] * m_factor
	t1['jackknifStackErrors'] = t1['jackknifStackErrors'] * m_factor
	t1['jackknifeSpectra'] = t1['jackknifeSpectra'] * m_factor
	t1['NspectraPerPixelMax'] = n.max(t1['NspectraPerPixel']) * n.ones_like(t1['NspectraPerPixel'])
	ts.append(t1)

t_out = Table( n.hstack(ts) )
t_out.remove_column('meanStack')
t_out.write(file_out, overwrite = True)

###########
###########
# T1 AGN ROSAT
###########
###########

T1_stacks = n.array( glob.glob( join(stack_dir,"X_AGN", "ROSAT_AGNT1*stack") ) )
T1_stacks.sort()
file_out =  join(stack_dir,"X_AGN", "ROSAT_AGNT1-stitched-stack.fits")

frac_Nspec = 0.5

idx_list = [0, 1, 11, 21, 31, 39]
wl_mins = [6900, 3200, 2500, 1800, 1480, 1150]
wl_maxs = [8100, 6900, 3200, 2500, 1800, 1480]

ii=-1
jj = idx_list[ii]
wl0, wl1  = wl_mins[ii], wl_maxs[ii]
print(T1_stacks[jj], wl0, wl1)
hdu = fits.open(T1_stacks[jj])[1].data
s_wl = (hdu['wavelength'] > wl0) & (hdu['wavelength'] < wl1) & (hdu['NspectraPerPixel'] > n.max(hdu['NspectraPerPixel']) * frac_Nspec)
med_norm = n.median(hdu['medianStack'][s_wl][-10:])
t0 = Table(hdu[s_wl])
t0['medianStack'] = t0['medianStack']/med_norm
t0['meanStack'] = t0['meanStack']/med_norm
t0['jackknifStackErrors'] = t0['jackknifStackErrors']/med_norm
t0['jackknifeSpectra'] = t0['jackknifeSpectra']/med_norm
t0['NspectraPerPixelMax'] = n.max(t0['NspectraPerPixel'])* n.ones_like(t0['NspectraPerPixel'])

ts = [ t0 ]
for jj, wl0, wl1 in zip(idx_list[::-1][1:], wl_mins[::-1][1:], wl_maxs[::-1][1:]):
	print(T1_stacks[jj], wl0, wl1)
	hdu = fits.open(T1_stacks[jj])[1].data
	s_wl = (hdu['wavelength']>wl0 ) & ( hdu['wavelength'] < wl1) & (hdu['NspectraPerPixel']>n.max(hdu['NspectraPerPixel']) * frac_Nspec )
	med_norm = n.median( hdu['medianStack'][s_wl][:10] )
	previous_norm = n.median( ts[-1]['medianStack'][-10:] )
	m_factor = previous_norm/med_norm
	t1 = Table(hdu[s_wl])
	t1['medianStack'] = t1['medianStack'] * m_factor
	t1['meanStack'] = t1['meanStack'] * m_factor
	t1['jackknifStackErrors'] = t1['jackknifStackErrors'] * m_factor
	t1['jackknifeSpectra'] = t1['jackknifeSpectra'] * m_factor
	t1['NspectraPerPixelMax'] = n.max(t1['NspectraPerPixel']) * n.ones_like(t1['NspectraPerPixel'])
	ts.append(t1)

t_out = Table( n.hstack(ts) )
t_out.remove_column('meanStack')
t_out.write(file_out, overwrite = True)



###########
###########
# T2 AGN ROSAT
###########
###########




T2_stacks = n.array( glob.glob( join(stack_dir,"X_AGN", "ROSAT_AGNT2*stack") ) )
T2_stacks.sort()
file_out =  join(stack_dir,"X_AGN", "ROSAT_AGNT2-stitched-stack.fits")

frac_Nspec = 0.5

idx_list = [0, 1, 7, 11]
wl_mins = [6100, 3650, 3000, 2300]
wl_maxs = [8100, 6100, 3650, 3000]

ii=-1
jj = idx_list[ii]
wl0, wl1  = wl_mins[ii], wl_maxs[ii]
hdu = fits.open(T2_stacks[jj])[1].data
s_wl = (hdu['wavelength'] > wl0) & (hdu['wavelength'] < wl1) & (hdu['NspectraPerPixel'] > n.max(hdu['NspectraPerPixel']) * frac_Nspec)
med_norm = n.median(hdu['medianStack'][s_wl][-10:])
t0 = Table(hdu[s_wl])
t0['medianStack'] = t0['medianStack']/med_norm
t0['meanStack'] = t0['meanStack']/med_norm
t0['jackknifStackErrors'] = t0['jackknifStackErrors']/med_norm
t0['jackknifeSpectra'] = t0['jackknifeSpectra']/med_norm
t0['NspectraPerPixelMax'] = n.max(t0['NspectraPerPixel'])* n.ones_like(t0['NspectraPerPixel'])

ts = [ t0 ]
for jj, wl0, wl1 in zip(idx_list[::-1][1:], wl_mins[::-1][1:], wl_maxs[::-1][1:]):
	hdu = fits.open(T2_stacks[jj])[1].data
	s_wl = (hdu['wavelength']>wl0 ) & ( hdu['wavelength'] < wl1) & (hdu['NspectraPerPixel']>n.max(hdu['NspectraPerPixel']) * frac_Nspec )
	med_norm = n.median( hdu['medianStack'][s_wl][:10] )
	previous_norm = n.median( ts[-1]['medianStack'][-10:] )
	m_factor = previous_norm/med_norm
	t1 = Table(hdu[s_wl])
	t1['medianStack'] = t1['medianStack'] * m_factor
	t1['meanStack'] = t1['meanStack'] * m_factor
	t1['jackknifStackErrors'] = t1['jackknifStackErrors'] * m_factor
	t1['jackknifeSpectra'] = t1['jackknifeSpectra'] * m_factor
	t1['NspectraPerPixelMax'] = n.max(t1['NspectraPerPixel']) * n.ones_like(t1['NspectraPerPixel'])
	ts.append(t1)

t_out = Table( n.hstack(ts) )
t_out.remove_column('meanStack')
t_out.write(file_out, overwrite = True)




list_2_stack = n.array(glob.glob(join(stack_dir, "*.ascii")))
for el in list_2_stack:
	stack_it(el)

for el in list_2_stack[::-1]:
	stack_it(el)

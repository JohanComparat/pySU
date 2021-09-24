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

from astropy.table import Table, Column

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115)#, Ob0=0.048206)
# create all input files :
path_2_elg3 = join(os.environ['HOME'],"SDSS/lss/catalogs/3", "inputs/ELG.v5_10_10.all.fits")
path_2_elg4 = join(os.environ['HOME'],"SDSS/lss/catalogs/4", "inputs/ELG.v5_11_0.rrv2.all.fits")
dr16_dir = join(os.environ['HOME'],"SDSS/dr16")
dr12_dir = join(os.environ['HOME'],"SDSS/dr12")

stack_dir = join(os.environ['HOME'],"SDSS/stacks")
p2_specObj = join(dr16_dir,"spAll-v5_13_0.fits")
p2_qso_3 = join(dr16_dir,"DR16Q_Superset_v3.fits")
p2_qso_4 = join(dr16_dir,"DR16Q_v4.fits")
p2_lrg = join(dr16_dir,"eBOSS_LRG_full_ALLdata-vDR16.fits")
p2_cmas = join(dr12_dir,"spAll-DR12.fits")

#p2_lrg_n = join(os.environ['HOME'],"data1/SDSS/eBOSS/lsscat/eBOSS_LRGpCMASS_clustering_data-NGC-vDR16.fits")
#p2_lrg_s = join(os.environ['HOME'],"data1/SDSS/eBOSS/lsscat/eBOSS_LRGpCMASS_clustering_data-SGC-vDR16.fits")
#p2_elg21 = join(os.environ['HOME'],"SDSS/catalogs/ELG_Y1.eboss21.fits")
#p2_elg22 = join(os.environ['HOME'],"SDSS/catalogs/ELG_Y1.eboss22.fits")
#p2_elg23 = join(os.environ['HOME'],"SDSS/catalogs/ELG_Y1.eboss23.fits")
#p2_qso = join(os.environ['HOME'],"data1/SDSS/eBOSS/lsscat/eBOSS_QSO_full_ALLdata-vDR16.fits")

spAll = fits.open(p2_specObj)[1].data
qso3 = fits.open(p2_qso_3)[1].data
#qso4 = fits.open(p2_qso_4)[1].data
#elg3 = fits.open(path_2_elg3)[1].data
elg4 = fits.open(path_2_elg4)[1].data
## QSO STACKS LISTS

cat = qso3[(qso3['Z']>0.05)&(qso3['Z']<6)]

Ngal = len(cat)
N_in_stack = 50000
N_factor = 1

#bins_2nd = n.arange(N_in_stack, N_in_stack*N_factor, N_in_stack)
print(Ngal)
#print(bins_2nd)

NNN,BBB=n.histogram(cat['Z'], bins=n.arange(0,7,0.001))
N_CM = n.cumsum(NNN)

N_bins = n.hstack(( n.arange(1, N_CM.max(), N_in_stack*N_factor), N_CM.max()-10 ))

itp = interp1d(N_CM, BBB[:-1]) 

z_mins = itp(N_bins)[:-1]
z_maxs = itp(N_bins)[1:]

def write_stack_list(cat, z_mins, z_maxs, prefix):
	plate = cat['PLATE']
	mjd = cat['MJD']              
	fiberid = cat['FIBERID']      
	z = cat['Z']
	for zmin, zmax in zip(z_mins, z_maxs):
		selection = (z > zmin) & (z < zmax)
		N_obj = len(selection.nonzero()[0])
		print(zmin, '<z<', zmax, 'N', N_obj)
		if N_obj>10:
			name = prefix+"_zmin_"+str(int(10*zmin)).zfill(2)+"_zmax_"+str(int(10*zmax)).zfill(2)+'.ascii'
			print(name)
			DATA = n.transpose([ plate, mjd, fiberid, z ])[selection]
			n.savetxt(os.path.join(stack_dir, name), DATA)

write_stack_list(cat, z_mins, z_maxs, 'STACK_DR16_QSO')

## ELG STACKS LISTS

cat = elg4[(elg4['rr_Z']>0.05) & (elg4['rr_Z']<1.7)] #

Ngal = len(cat)
N_in_stack = 30000
N_factor = 1

#bins_2nd = n.arange(N_in_stack, N_in_stack*N_factor, N_in_stack)
print(Ngal)
#print(bins_2nd)

NNN,BBB=n.histogram(cat['rr_Z'], bins=n.arange(0,7,0.001))
N_CM = n.cumsum(NNN)

N_bins = n.hstack(( n.arange(1, N_CM.max(), N_in_stack*N_factor), N_CM.max()-10 ))

itp = interp1d(N_CM, BBB[:-1])

z_mins = itp(N_bins)[:-1]
z_maxs = itp(N_bins)[1:]

def write_stack_list_elg(cat, z_mins, z_maxs, prefix):
	plate = cat['plate']
	mjd = cat['MJD']
	fiberid = cat['FIBERID']
	z = cat['rr_Z']
	for zmin, zmax in zip(z_mins, z_maxs):
		selection = (z > zmin) & (z < zmax) & (cat['rr_ZWARN']<=4)
		N_obj = len(selection.nonzero()[0])
		print(zmin, '<z<', zmax, 'N', N_obj)
		if N_obj>10:
			name = prefix+"_zmin_"+str(int(100*zmin)).zfill(2)+"_zmax_"+str(int(100*zmax)).zfill(2)+'.ascii'
			print(name)
			DATA = n.transpose([ plate, mjd, fiberid, z ])[selection]
			n.savetxt(os.path.join(stack_dir, name), DATA)

write_stack_list_elg(cat, z_mins, z_maxs, 'STACK_DR16_ELG')

## LRG STACK LISTS

lrg = fits.open(p2_lrg)[1].data
cma = fits.open(p2_cmas)[1].data

RG = ( (cma['BOSS_TARGET1'] & (1 | 2 | 3 | 7 | 8 )) > 0 )
cat1 = cma[(RG) & (cma['Z']>0.05) & (cma['Z']>cma['Z_ERR']) & (cma['ZWARNING']==0)]
cat2 = lrg [lrg['Z']>0.05]

cat = Table()
cat['PLATE'] = n.hstack(( cat1['PLATE'], cat2['PLATE'] ))
cat['MJD'] = n.hstack(( cat1['MJD'], cat2['MJD'] ))
cat['FIBERID'] = n.hstack(( cat1['FIBERID'], cat2['FIBERID'] ))
cat['Z'] = n.hstack(( cat1['Z'], cat2['Z'] ))

Ngal = len(cat)
N_in_stack = 100000
N_factor = 1

#bins_2nd = n.arange(N_in_stack, N_in_stack*N_factor, N_in_stack)
print(Ngal)
#print(bins_2nd)

NNN,BBB=n.histogram(cat['Z'], bins=n.arange(0.04,1.7,0.001))
N_CM = n.cumsum(NNN)

N_bins = n.arange(N_in_stack*N_factor, N_CM.max(), N_in_stack*N_factor)

itp = interp1d(N_CM, BBB[:-1])

z_mins = itp(N_bins)[:-1]
z_maxs = itp(N_bins)[1:]

def write_stack_list(cat, z_mins, z_maxs, prefix):
	plate = cat['PLATE']
	mjd = cat['MJD']
	fiberid = cat['FIBERID']
	z = cat['Z']
	for zmin, zmax in zip(z_mins, z_maxs):
		selection = (z > zmin) & (z < zmax)
		N_obj = len(selection.nonzero()[0])
		print(zmin, '<z<', zmax, 'N', N_obj)
		if N_obj>10:
			name = prefix+"_zmin_"+str(int(100*zmin)).zfill(3)+"_zmax_"+str(int(100*zmax)).zfill(3)+'.ascii'
			print(name)
			DATA = n.transpose([ plate, mjd, fiberid, z ])[selection]
			n.savetxt(os.path.join(stack_dir, name), DATA)

write_stack_list(cat, z_mins, z_maxs, 'STACK_DR16_LRG')

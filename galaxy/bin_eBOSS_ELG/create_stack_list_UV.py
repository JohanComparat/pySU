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


from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115)#, Ob0=0.048206)
# create all input files :
path_2_cat = join(os.environ['HOME'],"SDSS/lss/catalogs/3", "inputs/ELG.v5_10_10.all.fits")

#path_2_cat = join(os.environ['HOME'],"SDSS/lss/catalogs/4", "inputs/ELG.v5_11_0.rrv2.all.fits")

cat = fits.open(path_2_cat)[1].data

Ngal = len(cat)
N_in_stack = 10000
N_factor = 4

#bins_2nd = n.arange(N_in_stack, N_in_stack*N_factor, N_in_stack)
print(Ngal)
#print(bins_2nd)

NNN,BBB=n.histogram(cat['Z'], bins=n.arange(0,4,0.001))
N_CM = n.cumsum(NNN)

N_bins = n.arange(N_in_stack*N_factor, N_CM.max(), N_in_stack*N_factor)

itp = interp1d(N_CM, BBB[:-1]) 

z_mins = itp(N_bins)[:-1]
z_maxs = itp(N_bins)[1:]


# CREATES A few stacks as a function of [OII] EW

z0,z1 = 0.6, 1.2
z_sel = (cat['rr_Z']>z0) & (cat['rr_Z']<z1) & (cat['rr_ZWARN']<=4)

qty = 'rr_Z_O2_3728_EW'
qtyN = 'O2EW'
if qtyN == 'O2EW':
  sel = (z_sel)&(cat[qty]>0.1)&(cat[qty]<10000)

NNN,BBB=n.histogram(n.log10(cat[qty][sel]), bins=100000)
N_CM = n.cumsum(NNN)
N_bins = n.arange(N_in_stack*N_factor, N_CM.max(), N_in_stack*N_factor)
itp = interp1d(N_CM, 10**BBB[:-1]) 
qty_mins = itp(N_bins)[:-1]
qty_maxs = itp(N_bins)[1:]
for q0,q1 in zip(qty_mins, qty_maxs):
	print(q0,q1)
	selection = ( sel ) & (cat[qty]>=q0)&(cat[qty]<q1)
	DATA = n.transpose([ cat['plate'], cat['MJD'], cat['FIBERID'], cat['rr_Z'] ]) [selection]
	path_2_input = join(os.environ['HOME'],"SDSS/stacks", "eboss-elg_"+str(z0)+"_z_"+str(z1)+"_"+str(q0)+"_"+qtyN+"_"+str(q1)+".asc")
	print(path_2_input)
	print(len(DATA))
	n.savetxt(path_2_input, DATA)

# CREATES A few stacks as a function of [OII] L
sel = (z_sel) & (cat['rr_Z_O2_3728_flux']>0) & (cat['rr_Z_O2_3728_flux']<1)
dL = cosmoMD.luminosity_distance(cat['rr_Z'][sel]).to(u.cm).value
Loii = n.log10(cat['rr_Z_O2_3728_flux'][sel]*4*n.pi*dL**2)
qtyN = 'O2lum' 
NNN,BBB=n.histogram(Loii, bins=100000)
N_CM = n.cumsum(NNN)
N_bins = n.arange(N_in_stack*N_factor, N_CM.max(), N_in_stack*N_factor)
itp = interp1d(N_CM, BBB[:-1]) 
qty_mins = itp(N_bins)[:-1]
qty_maxs = itp(N_bins)[1:]
for q0,q1 in zip(qty_mins, qty_maxs):
	print(q0,q1)
	selection =  (Loii>=q0)&(Loii<q1)
	DATA = n.transpose([ cat['plate'][sel], cat['MJD'][sel], cat['FIBERID'][sel], cat['rr_Z'][sel] ]) [selection]
	path_2_input = join(os.environ['HOME'],"SDSS/stacks", "eboss-elg_"+str(z0)+"_z_"+str(z1)+"_"+str(q0)+"_"+qtyN+"_"+str(q1)+".asc")
	print(path_2_input)
	print(len(DATA))
	n.savetxt(path_2_input, DATA)


def create_lists(qty, qtyN):
	for z0,z1 in zip(z_mins, z_maxs):
		print('-------------------------------------')
		print(z0,z1)
		z_sel = (cat['rr_Z']>z0) & (cat['rr_Z']<z1)
		bins_2nd = n.array([ int(len(cat['rr_Z'][z_sel])/3.), int(len(cat['rr_Z'][z_sel])*2./3.) ])
		if qtyN == 'O2flux' :
			z_sel = (cat['rr_Z']>z0) & (cat['rr_Z']<z1) & (cat[qty]>0)
			bins_2nd = n.array([ int(len(cat['rr_Z'][z_sel])/3.), int(len(cat['rr_Z'][z_sel])*2./3.) ])
		if qtyN == 'rw1' :
			z_sel = (cat['rr_Z']>z0) & (cat['rr_Z']<z1) & (cat[qty]>-10)
			bins_2nd = n.array([ int(len(cat['rr_Z'][z_sel])/3.), int(len(cat['rr_Z'][z_sel])*2./3.) ])
		print(len(cat[qty][z_sel]))
		print(bins_2nd)
		qty_bins = n.arange(n.min(cat[qty][z_sel]), n.max(cat[qty][z_sel]), (-n.min(cat[qty][z_sel]) + n.max(cat[qty][z_sel]))/1000000. )
		itp2 = interp1d( n.cumsum(n.histogram(cat[qty][z_sel], bins=qty_bins)[0]), qty_bins[:-1] )
		print(itp2.x,itp2.y)
		qty_mins = n.hstack((n.min(cat[qty][z_sel]), itp2(bins_2nd) ))
		qty_maxs =  n.hstack(( itp2(bins_2nd), n.max(cat[qty][z_sel]) ))
		for q0,q1 in zip(qty_mins, qty_maxs):
			print(q0,q1)
			selection = ( z_sel ) & (cat[qty]>q0)&(cat[qty]<q1)
			DATA = n.transpose([ cat['plate'], cat['MJD'], cat['FIBERID'], cat['rr_Z'] ]) [selection]
			path_2_input = join(os.environ['HOME'],"SDSS/stacks", "eboss-elg_"+str(z0)+"_z_"+str(z1)+"_"+str(q0)+"_"+qtyN+"_"+str(q1)+".asc")
			print(path_2_input)
			print(len(DATA))
			n.savetxt(path_2_input, DATA)

qty = 'rr_Z_O2_3728_EW'
qtyN = 'O2EW'

#create_lists(qty = 'rr_Z_O2_3728_flux' , qtyN = 'O2flux'  )
create_lists(qty = 'rw1'               , qtyN = 'rw1'     )
create_lists(qty = 'rr_fast_lmass'     , qtyN = 'mass'    )
create_lists(qty = 'g'                 , qtyN = 'g'       )
create_lists(qty = 'gr'                , qtyN = 'gr'      )
create_lists(qty = 'rz'                , qtyN = 'rz'      )

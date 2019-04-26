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
#path_2_cat = join(os.environ['HOME'],"SDSS/lss/catalogs/3", "inputs/ELG.v5_10_10.all.fits")

path_2_cat = join(os.environ['HOME'],"SDSS/lss/catalogs/4", "inputs/ELG.v5_11_0.rrv2.all.fits")

cat = fits.open(path_2_cat)[1].data

Ngal = len(cat)
N_in_stack = 200000
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

z0,z1 = 0.2, 1.5
selection = (cat['rr_Z']>z0) & (cat['rr_Z']<z1) & (cat['rr_ZWARN']<=4)

ids_sort = n.argsort(cat['rr_Z'][selection])

DATA = n.transpose([ cat['plate'], cat['MJD'], cat['FIBERID'], cat['rr_Z'] ]) [selection][ids_sort]
path_2_input = join(os.environ['HOME'],"SDSS/stacks", "eboss-elg_"+str(z0)+"_z_"+str(z1)+".asc")
print(path_2_input)
print(len(DATA))
n.savetxt(path_2_input, DATA)

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
path_2_cat = os.path.join(os.environ['HOME'], 'wwwDir/stuff/catalogue_qso_full.fits')
#path_2_cat = join(os.environ['HOME'],"SDSS/lss/catalogs/3", "inputs/ELG.v5_10_10.all.fits")
'SDSS/stack_cluster_qso_pairs'
#path_2_cat = join(os.environ['HOME'],"SDSS/lss/catalogs/4", "inputs/ELG.v5_11_0.rrv2.all.fits")

cat_i = fits.open(path_2_cat)[1].data
selection = (cat_i['cluster_z']>0.3)&(cat_i['cluster_z']<0.4)&(cat_i['galaxy_z']>cat_i['cluster_z'])
cat = cat_i[selection]

Ngal = len(cat)
N_in_stack = 10000
N_factor = 4

f_radius = cat['angular_separation']/cat['cluster_r200c_deg']
z_corr = (1+cat['galaxy_z'])/(1+cat['cluster_z'])-1

b1 = (f_radius<=0.5)
b2 = (f_radius>0.5)&(f_radius<=1.)
b3 = (f_radius>1.0)&(f_radius<=2.0)
b4 = (f_radius>2.0)&(f_radius<=3.0)

DATA = n.transpose([ cat['plate'], cat['MJD'], cat['FIBERID'], z_corr ]) [b1]
path_2_input = join(os.environ['HOME'],"SDSS/stack_cluster_qso_pairs", "clusterXqso_b1.asc")
print(path_2_input, len(DATA))
n.savetxt(path_2_input, DATA)

DATA = n.transpose([ cat['plate'], cat['MJD'], cat['FIBERID'], z_corr ]) [b2]
path_2_input = join(os.environ['HOME'],"SDSS/stack_cluster_qso_pairs", "clusterXqso_b2.asc")
print(path_2_input, len(DATA))
n.savetxt(path_2_input, DATA)

DATA = n.transpose([ cat['plate'], cat['MJD'], cat['FIBERID'], z_corr ]) [b3]
path_2_input = join(os.environ['HOME'],"SDSS/stack_cluster_qso_pairs", "clusterXqso_b3.asc")
print(path_2_input, len(DATA))
n.savetxt(path_2_input, DATA)

DATA = n.transpose([ cat['plate'], cat['MJD'], cat['FIBERID'], z_corr ]) [b4]
path_2_input = join(os.environ['HOME'],"SDSS/stack_cluster_qso_pairs", "clusterXqso_b4.asc")
print(path_2_input, len(DATA))
n.savetxt(path_2_input, DATA)

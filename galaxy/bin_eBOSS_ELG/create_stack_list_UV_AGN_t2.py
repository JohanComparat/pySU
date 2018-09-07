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


# create all input files :

path_2_cat = os.path.join(os.environ['HOME'], 'SDSS/dr14/VAC_2RXS_DR14_ALLW_SDSS_GAIA_FIRST_final_SPM.FITS')

cat = fits.open(path_2_cat)[1].data

z_bin = (cat['VI_Z']>0.1)&(cat['VI_Z']<1.2)

selection = (z_bin)&((cat['VI_CLASS']=='NLAGN') | (cat['VI_CLASS']=='GAL'))

DATA = n.transpose([ cat['plate'], cat['MJD'], cat['FIBERID'], cat['Z'] ]) [selection]
path_2_input = join(os.environ['HOME'],"SDSS/lss/catalogs/3", "xagn_0.38.asc")
n.savetxt(path_2_input, DATA)

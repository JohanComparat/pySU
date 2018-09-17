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

z_bin = (cat['VI_Z']>0.0)#&(cat['VI_Z']<1.2)


selection = (z_bin)&((cat['VI_CLASS']=='NLAGN') | (cat['VI_CLASS']=='GAL'))

DATA0 = n.transpose([ cat['plate'], cat['MJD'], cat['FIBERID'], cat['Z'] ]) [selection]



path_2_s82 = os.path.join(os.environ['OBS_REPO'], 'stripe82X' )
path_2_cat = os.path.join(path_2_s82, 'final_catalog_updated_dr16_v51010.fits')
cat = fits.open(path_2_cat)[1].data

#
zspec= (cat['Z_v51010']>0)#&(cat['Z_v51010']<1.5)
# template selections
stars   = (cat['MORPHOLOGY'] == 1 ) & ( zspec) #& (cat['r']>0)
gal_ell = (cat['MORPHOLOGY'] == 2 ) & ( zspec) #& (cat['r']>0)
gal_spi = (cat['MORPHOLOGY'] == 3 ) & ( zspec) #& (cat['r']>0)
type2   = (cat['MORPHOLOGY'] == 4 ) & ( zspec) #& (cat['r']>0)
gal_SB  = (cat['MORPHOLOGY'] == 5 ) & ( zspec) #& (cat['r']>0)
type1   = (cat['MORPHOLOGY'] == 6 ) & ( zspec) #& (cat['r']>0)
qso     = ( (cat['MORPHOLOGY'] == 6 ) | (cat['MORPHOLOGY'] == 7 )) & ( zspec) #& (cat['r']>0)


selection = type2

DATA = n.transpose([ cat['plate_v51010'], cat['MJD_v51010'], cat['FIBERID_v51010'], cat['Z_v51010'] ]) [selection]

for el in DATA:
  isin = (DATA0.T[0]==el[0])&(DATA0.T[2]==el[2])&(DATA0.T[1]==el[1])
  print(el, DATA0[isin])
  if len(isin.nonzero()[0])==0:
    print('new')
    DATA0 = n.vstack((DATA0,el))

path_2_input = join(os.environ['HOME'],"SDSS/lss/catalogs/3", "s82xagn.asc")
n.savetxt(path_2_input, DATA)

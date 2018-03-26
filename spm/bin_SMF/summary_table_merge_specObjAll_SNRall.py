from astropy.io import fits
import numpy as n
import glob
import os
import sys

vv = sys.argv[1]

#plates = n.array(os.listdir("/home/comparat/SDSS/26/stellarpop"))
top_dir = "/home/comparat/SDSS/"+vv
snr_dir = top_dir+"/SNR"
cat_dir = top_dir+"/catalogs"

snr_file = snr_dir + "/snr-summary-wh.ascii"

if vv=='26':
  #cat_file = cat_dir + "/specObj-SDSS-dr14.fits"
  cat_file = cat_dir + "/spFlyAll_906.fits"
  #cat_file = cat_dir + "/spFlyAll_672.fits"

if vv=='v5_10_0':
  #cat_file = cat_dir + "/specObj-BOSS-dr14.fits"
  cat_file = cat_dir + "/spFlyAll_906.fits"
  

cm_3d = lambda f1, f2, out : """stilts tmatch2 in1="""+f1+""" ifmt1=fits in2="""+f2+""" ifmt2=ascii matcher=3d values1="PLATE MJD FIBERID" values2="PL MJ FI" params=0.0001 join=all1 find=best ocmd='sort plate' ocmd='delcols "PL MJ FI"' omode=out ofmt=fits out="""+out


out = cat_dir+'/FireFly_672.fits'
c1 = cm_3d(cat_file, snr_file, out)
print(c1)
os.system(c1)

# PL MJ FI SNR_ALL SNR_32_35 SNR_35_39 SNR_39_41 SNR_41_55 SNR_55_68 SNR_68_74 SNR_74_93


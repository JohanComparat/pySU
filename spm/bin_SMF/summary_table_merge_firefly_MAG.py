from astropy.io import fits
import numpy as n
import glob
import os
import sys

vv = sys.argv[1]

#plates = n.array(os.listdir("/home/comparat/SDSS/26/stellarpop"))
top_dir = "/data36s/comparat/SDSS/"+vv
cat_dir = top_dir+"/catalogs"
cat_file  = cat_dir+'/FireFly.fits'
#cat_file  = cat_dir+'/FireFly_672.fits'

mag_file = "/data36s/comparat/SDSS/specphotalldr14_comparat.fit"

cm_3d = lambda f1, f2, out : """stilts tmatch2 in1="""+f1+""" ifmt1=fits in2="""+f2+""" ifmt2=fits matcher=3d values1="PLATE MJD FIBERID" values2="PLATE MJD FIBERID" params=0.0001 join=all1 find=best omode=out ofmt=fits out="""+out

#out = cat_dir+'/FireFly_672_mag.fits'
out = cat_dir+'/FireFly_mag.fits'

c1 = cm_3d(cat_file, mag_file, out)
print(c1)
os.system(c1)


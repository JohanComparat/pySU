from astropy.io import fits
import numpy as n
import glob
import os
import sys

vv = sys.argv[1]

#plates = n.array(os.listdir("/home/comparat/SDSS/26/stellarpop"))
top_dir = "/home/comparat/SDSS/"+vv+"/stellarpop"

command = lambda lis, outname : """stilts tcat in=@"""+lis+""" ifmt=fits ocmd='sort plate' omode=out ofmt=fits out="""+outname

lis = top_dir+'/flyAllList_906'
out = "/home/comparat/SDSS/"+vv+ "/catalogs/spFlyAll_906.fits"
c1 = command(lis, out)
print(c1)
os.system(c1)

lis = top_dir+'/flyAllList_159'
out = "/home/comparat/SDSS/"+vv+ "/catalogs/spFlyAll_159.fits"
c1 = command(lis, out)
print(c1)
os.system(c1)

import glob
from os.path import join
import numpy as n
import astropy.io.fits as fits
import lib_functions_1pt as lib
import os

#Quantity studied
qty = "mvir"

# General information
zList_all =  join(os.environ['PYSU_MD_DIR'], "data", "z-list-all-boxes.txt") 
z0 = n.loadtxt(zList_all,unpack=True)
zList_all2 =  join(os.environ['PYSU_MD_DIR'], "data", "z-list-2LINEAR-COSMO.txt") 
z0short = n.loadtxt(zList_all2,unpack=True,dtype='S')

# redshift lists
dir_boxes =  n.array([os.environ['MD04_DIR'], os.environ['MD10_DIR'], os.environ['MD25_DIR'], os.environ['MD40_DIR'], os.environ['MD25NW_DIR'], os.environ['MD40NW_DIR']])
zList_files = n.array([ join(dir_box,"redshift-list.txt") for dir_box in dir_boxes])

# one point function lists
fileC = n.array(glob.glob( join(os.environ['MULTIDARK_LIGHTCONE_DIR'],"MD_*Gpc*", "properties", qty,"*t_*_Central_JKresampling.pkl")))
fileB = n.array(glob.glob( join( os.environ['MULTIDARK_LIGHTCONE_DIR'],"MD_*Gpc*","properties", qty,"*t_*_"+qty+"_JKresampling.bins")))
fileS = n.array(glob.glob( join( os.environ['MULTIDARK_LIGHTCONE_DIR'],"MD_*Gpc*","properties", qty,"*t_*_Satellite_JKresampling.pkl")))

print "considers ",len(fileC), qty , " function files"


for ii, el in enumerate(fileC):
	print el
	print fileS[ii]
	print fileB[ii]
	lib.convert_pkl_mass(fileC[ii], fileS[ii], fileB[ii],zList_files, z0, z0short, qty)

af = n.array(glob.glob(join(os.environ['MULTIDARK_LIGHTCONE_DIR'], qty, "data", "MD_*_"+qty+".fits") ) )
d0 = fits.open(af[0])[1].data
for ii in range(1,len(af),1):
	d1 = fits.open(af[ii])[1].data
	d0 = n.hstack((d0,d1))

hdu2 = fits.BinTableHDU.from_columns(d0)
os.system("rm "+join(os.environ['MULTIDARK_LIGHTCONE_DIR'], qty, "MD_"+qty+"_summary.fits"))
hdu2.writeto( join(os.environ['MULTIDARK_LIGHTCONE_DIR'], qty, "MD_"+qty+"_summary.fits") )




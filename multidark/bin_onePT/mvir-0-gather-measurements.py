import glob
from os.path import join
import numpy as n
import astropy.io.fits as fits
import lib_functions_1pt as lib
import os


#Quantity studied
qty = "mvir"

# redshift lists
dir_boxes =  n.array([os.environ['MD04_DIR'], os.environ['MD10_DIR'], os.environ['MD25_DIR'], os.environ['MD40_DIR'], os.environ['MD25NW_DIR'], os.environ['MD40NW_DIR']])
zList_files = n.array([ join(dir_box,"redshift-list.txt") for dir_box in dir_boxes])

# one point function lists
fileC = n.array(glob.glob( join(os.environ['MULTIDARK_LIGHTCONE_DIR'],"MD_*Gpc*", "properties", qty,"*t_*_Central_JKresampling.pkl")))
fileB = n.array(glob.glob( join( os.environ['MULTIDARK_LIGHTCONE_DIR'],"MD_*Gpc*","properties", qty,"*t_*_"+qty+"_JKresampling.bins")))
fileS = n.array(glob.glob( join( os.environ['MULTIDARK_LIGHTCONE_DIR'],"MD_*Gpc*","properties", qty,"*t_*_Satellite_JKresampling.pkl")))
"""
fileC = n.array(glob.glob( join(os.environ['MULTIDARK_LIGHTCONE_DIR'],"MD_4GpcNW", "properties", qty,"*t_*_Central_JKresampling.pkl")))
fileB = n.array(glob.glob( join( os.environ['MULTIDARK_LIGHTCONE_DIR'],"MD_4GpcNW","properties", qty,"*t_*_"+qty+"_JKresampling.bins")))
fileS = n.array(glob.glob( join( os.environ['MULTIDARK_LIGHTCONE_DIR'],"MD_4GpcNW","properties", qty,"*t_*_Satellite_JKresampling.pkl")))
"""
fileC.sort()
fileS.sort()
fileB.sort()
print len(fileC), len(fileB), len(fileS)

print "considers ",len(fileC), qty , " function files"

for ii, el in enumerate(fileC):
	print el
	print fileS[ii]
	print fileB[ii]
	lib.convert_pkl_mass(fileC[ii], fileS[ii], fileB[ii],zList_files, qty)

af = n.array(glob.glob(join(os.environ['MULTIDARK_LIGHTCONE_DIR'], qty, "data", "MD_*_"+qty+".fits") ) )
d0 = fits.open(af[0])[1].data
for ii in range(1,len(af),1):
	d1 = fits.open(af[ii])[1].data
	d0 = n.hstack((d0,d1))

hdu2 = fits.BinTableHDU.from_columns(d0)

writeName = join(os.environ['MULTIDARK_LIGHTCONE_DIR'], qty, "MD_"+qty+"_summary.fits")
if os.path.isfile(writeName):
	os.remove(writeName)
	
hdu2.writeto( writeName )




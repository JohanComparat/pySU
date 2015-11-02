import time
t0 = time.time()
import numpy as n
import os
from os.path import join
from HaloSelection import *

zmin_IN, zmax_IN, nGal_Deg2_IN  = n.loadtxt("/data2/DATA/eBOSS/ELG/observations/NZ/nz-fisherGRIW1.dat",unpack=True)
ok = (zmin_IN >= 0.4) & (zmax_IN <= 1.2) 
zmin, zmax, nGal_Deg2 = zmin_IN[ok], zmax_IN[ok], nGal_Deg2_IN[ok] 

mockOutput_dir = "/data2/DATA/eBOSS/ELG-fischer-mocks/"
#"/Volumes/data/MD-lightcones/mocks-v1/"

lcDir = "/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/lightcones/lc_square_0.1z1.4/"
lcFile = join(lcDir,"lightcone_MD_2.5Gpc_0.4z1.4.fits")
area = (2*30.)**2.


hdu = fits.open(lcFile)
mockName="xx"

mm = MultiDarkMock(hdu, area, mockOutput_dir, mockName, zmin, zmax, nGal_Deg2)

mm.initialize()
t1 = time.time()
print "dt=",t1 - t0

mm.mockName = "tryMocks-sham"
mm.make_sham_catalog(colN='mvir')
mm.write_full_catalog_fits()
mm.write_catalog_ascii()

t1 = time.time()
mm.get_distrib_QTY('mvir',0.5,0.6)
t2 = time.time()
print "dt=",t2 - t1

# 30 + 180

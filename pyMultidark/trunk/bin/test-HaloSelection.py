import time
t0 = time.time()
import numpy as n
import os
from os.path import join
from HaloSelection import *

zmin, zmax, nGal_Deg2  =  n.arange(0.4,1.2,0.1), n.arange(0.5,1.31,0.1), n.ones_like(n.arange(0.4,1.2,0.1))*200
mockOutput_dir = "/Volumes/data/MD-lightcones/mocks-v1/"
lcDir = mockOutput_dir
lcFile = join(lcDir,"lightcone_MD_2.5Gpc_0.4z1.4.fits")
area = (2*30.)**2.


hdu = fits.open("/Volumes/data/MD-lightcones/mocks-v1/lightcone_MD_2.5Gpc_0.4z1.4.fits")
mockName="xx"

mm = MultiDarkMock(hdu, area, mockOutput_dir, mockName, zmin, zmax, nGal_Deg2)

mm.initialize()
t1 = time.time()
print t1 - t0

mm.mockName = "tryMocks-sham"
mm.make_sham_catalog(colN='mvir')
mm.write_full_catalog_fits()
mm.write_catalog_ascii()
mm.get_distrib_QTY('mvir',0.5,0.6)
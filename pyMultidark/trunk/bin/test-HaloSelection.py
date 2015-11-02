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

t2 = time.time()
print "dt=",t2 - t1
t3 = time.time()
print "dt=",t3 - t2

sats=n.arange(0.,0.31,0.025)
mMean=n.arange(12.2,13.6,0.2)
for ii in range(len(sats)):
    for jj in range(len(mMean)):
        p1 = mMean[jj]
        p2 = mMean[jj]-0.2
        p3 = sats[ii]
        mm.mockName = "tryMocks-gaussian_mean_"+str(p1)+"_sig_"+str(p2)+"_fsat_"+str(p3)
        mm.make_GaussianFsat_catalog("mvir", n.ones_like(zmin)*10**p1,n.ones_like(zmin)*10**p2, n.ones_like(zmin)*p3)
        #mm.write_full_catalog_fits()
        mm.write_catalog_ascii()
        mm.create_random_catalog()
        mm.writeClusteringParamFile("monopole")
        #mm.writeClusteringParamFile("angular","_d1")
        mm.writeClusteringParamFile("angular","_d2")
        mm.writeClusteringParamFile("angular","_d3")
        mm.compute_clustering()

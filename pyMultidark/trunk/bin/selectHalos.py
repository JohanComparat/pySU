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


sats=n.arange(0.,0.31,0.025)
mMean=n.arange(12.2,13.6,0.2)
"""
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

"""

data_clustering_dir= '/data2/DATA/eBOSS/ELG/observations/clustering-fischer/'
aa,bb,cc=n.loadtxt(join(data_clustering_dir,'grizw.xi2'),unpack=True,usecols=(0,1,2))

x_data = [aa,bb,2.*cc]

xx0,yy0,y0E=n.loadtxt(join(data_clustering_dir,'elg-fischer-griw-all-2PCF-angular-d3.dat'),unpack=True,usecols=(0,1,2))
xx1,yy1,y1E=n.loadtxt(join(data_clustering_dir,'elg-fischer-griw-all-2PCF-angular-d2.dat'),unpack=True,usecols=(0,1,2))
xx2,yy2,y2E=n.loadtxt(join(data_clustering_dir,'elg-fischer-griw-all-2PCF-angular-d1.dat'),unpack=True,usecols=(0,1,2))
xx=n.hstack((xx0,xx1[1:],xx2[1:]))
yy=n.hstack((yy0,yy1[1:],yy2[1:]))
yE=n.hstack((y0E,y1E[1:],y2E[1:]))*2.

w_data = [xx,yy,yE]
chi2map=[]
for ii in range(len(sats)):
    for jj in range(len(mMean)):
        p1 = mMean[jj]
        p2 = mMean[jj]-0.2
        p3 = sats[ii]
        mm.mockName = "tryMocks-gaussian_mean_"+str(p1)+"_sig_"+str(p2)+"_fsat_"+str(p3)
        chi2Wr, chi2Xr = mm.compare_clustering_data_mock(w_data, x_data, s_max_chi2=19)
        chi2map.append([p1,p2,p3,chi2Wr, chi2Xr])
        

chi2map = n.transpose(chi2map)

p.figure(10,(6,5))
p.axes([0.18,0.15,0.75,0.75])
p.ylabel('f sat')
p.xlabel('M mean')
#p.plot( chi2map[3], chi2map[1][chiC<1.5], 'k*',mec='w',ms=20,alpha=0.2)
p.scatter(chi2map[2], chi2map[0],c=(chi2map[3]+chi2map[4])/2.,s=60)
cb=p.colorbar()
cb.set_label('chi2 / ndof')
p.savefig(join(data_clustering_dir,"Fischer-griw-Mmean-fsat.pdf"))
p.clf()
    
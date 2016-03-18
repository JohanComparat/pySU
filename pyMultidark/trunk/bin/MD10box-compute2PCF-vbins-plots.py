import astropy.units as uu
import astropy.cosmology as co
aa = co.Planck13
#aah = co.FlatLambdaCDM(H0=100.0 *uu.km / (uu.Mpc *uu.s), Om0=0.307, Tcmb0=2.725 *uu.K, Neff=3.05, m_nu=[ 0.  ,  0. ,   0.06]*uu.eV, Ob0=0.0483)
#rhom = aa.critical_density0.to(uu.solMass*uu.Mpc**-3).value
#aa.critical_density0.to(uu.solMass*uu.Mpc**-3).value
#aah.critical_density0.to(uu.solMass*uu.Mpc**-3).value

from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp1d
import numpy as n
import matplotlib
#matplotlib.use('pdf')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p
import glob
import sys
from scipy.optimize import curve_fit
import cPickle
from os.path import join
from scipy.optimize import minimize

dir = join("D:","data","MultiDark")

dir_04 = join(dir,"MD_0.4Gpc")
dir_10 = join(dir,"MD_1Gpc")
dir_25 = join(dir,"MD_2.5Gpc")
dir_40 = join(dir,"MD_4Gpc")

dir_04_b = join(dir,"MD_0.4Gpc","halo_bias","clustering")
dir_10_b = join(dir,"MD_1Gpc","halo_bias","clustering")
dir_25_b = join(dir,"MD_2.5Gpc","halo_bias","clustering")
dir_40_b = join(dir,"MD_4Gpc","halo_bias","clustering")

dir_boxes =  n.array([dir_04, dir_10, dir_25, dir_40])
zList_files = n.array([ join(dir_box,"redshift-list.txt") for dir_box in dir_boxes])

snn, scaleFactor = n.loadtxt(join(dir_10,"redshift-list.txt"),unpack=True)
#zList_files[2],unpack=True)

sn2a = {scaleFactor[i]: redshift[i] for i in range(len(snn))}

ll = n.array( glob.glob(join(dir_10,"halo_bias","clustering","*.pkl")))#dir_25_b,"*.pkl" ) ))
ll.sort()
scaleFact = n.empty(len(ll))
vmin = n.empty(len(ll))
vmax = n.empty(len(ll))
zz = n.empty(len(ll))
for ii in range(len(ll)):
	scaleFact[ii] = float(ll[ii].split('\\')[-1].split('_')[1])
	vmin[ii] = float(ll[ii].split('\\')[-1].split('_')[3])
	vmax[ii] = float(ll[ii].split('\\')[-1].split('_')[4])
	zz[ii] = 1./scaleFact[ii] -1.

f=open("D:\data\MultiDark\MD_2.5Gpc\halo_bias\clustering\hlist_46_vmax_200.0_224.4_xiR.pkl",'r')
bin_xi3D_25, xi_25_200 = cPickle.load(f)
f.close()
f=open("D:\data\MultiDark\MD_2.5Gpc\halo_bias\clustering\hlist_46_vmax_224.4_251.79_xiR.pkl",'r')
bin_xi3D_25, xi_25_224 = cPickle.load(f)
f.close()
rr_25 = (bin_xi3D_25[1:] + bin_xi3D_25[:-1])/2.


sel = (zz>2.9)
xi = n.empty((len(ll[sel]),20))
for ii, file in enumerate(ll[sel]):
	f=open(file,'r')
	bin_xi3D, xi[ii] = cPickle.load(f)
	f.close()

rr = (bin_xi3D[1:] + bin_xi3D[:-1])/2.
	
p.figure(0,(10,6))
p.axes([0.17,0.17,0.6,0.8])
p.plot(rr_25,xi_25_200,ls='dashed', lw=2,label="L2.5 z0.3 vmax200 ")
#p.plot(rr_25,xi_25_224,ls='dashed', lw=2,label="L2.5 z0.3 vmax224 ")
for ii in n.arange(len(ll[sel])):
	p.plot(rr,xi[ii],label=str(n.round(vmin[sel][ii])))

p.xlabel('r Mpc/h')
p.ylabel(' r xi(r)')
p.yscale('log')
p.xscale('log')
gl = p.legend(bbox_to_anchor=(1.05, 1), loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.show()
	
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

snn, redshift = n.loadtxt("redshift-list.txt",unpack=True)
#zList_files[2],unpack=True)

sn2z = {snn[i]: redshift[i] for i in range(len(snn))}

ll = n.array( glob.glob(join("halo_bias","clustering","*.pkl")))#dir_25_b,"*.pkl" ) ))
ll.sort()
snapNum = n.empty(len(ll),dtype=int)
vmin = n.empty(len(ll))
zz = n.empty(len(ll))
for ii in range(len(ll)):
	snapNum[ii] = int(ll[ii].split('\\')[-1].split('_')[1])
	vmin[ii] = float(ll[ii].split('\\')[-1].split('_')[4])
	zz[ii] = sn2z[snapNum[ii]]


sel = (snapNum ==46)
xi = n.empty((len(ll[sel]),19))
for ii, file in enumerate(ll[sel]):
	f=open(file,'r')
	bin_xi3D, xi[ii] = cPickle.load(f)
	f.close()

rr = (bin_xi3D[1:] + bin_xi3D[:-1])/2.
	
p.figure(0,(10,6))
p.axes([0.17,0.17,0.6,0.8])
for ii in n.arange(len(ll[sel]))[::3]:
	p.plot(rr,rr*rr*xi[ii],label=str(n.round(vmin[sel][ii])))

p.xlabel('r Mpc/h')
p.ylabel(' r xi(r)')
#p.yscale('log')
#p.xscale('log')
gl = p.legend(bbox_to_anchor=(1.05, 1), loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.show()
	
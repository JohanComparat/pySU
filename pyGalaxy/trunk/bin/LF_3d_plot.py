import numpy as n
import matplotlib
matplotlib.rcParams['font.size']=14
import matplotlib.pyplot as p
import glob
import sys
from scipy.optimize import curve_fit
import cPickle
import os
from os.path import join
data_dir = os.environ['DATA_DIR']

Pdir = join(data_dir,"Products_Galaxies", "emissionLineLuminosityFunctions")
# "/home/comparat/database/Products_Galaxies/emissionLineLuminosityFunctions/" # on eboss
lines = "H1_4862", "O3_5007", "O2_3728"

###################### O2 3728 #######################
###################### O2 3728 #######################
###################### O2 3728 #######################
line = lines[2]
files = n.array(glob.glob(join(Pdir,line, "*.points")))
files.sort()
print files, line

logls = lambda z : -0.182 * z + 41.272
ps = lambda z: 0.741 * z - 2.78
a = -1.8
sig = 0.53
saundersFct=lambda logl, z : 10**ps(z) * (10**logl/10**logls(z))**(a+1) * n.e**( -n.log10( 1 +10**logl/10**logls(z))**2./(2*sig**2.))

loglsUP = lambda z : -0.182 * z + 41.421
psUP = lambda z: 0.741 * z - 3.026
saundersFctUp=lambda logl, z : 10**psUP(z) * (10**logl/10**loglsUP(z))**(a+1) * n.e**( -n.log10( 1 +10**logl/10**loglsUP(z))**2./(2*sig**2.))

loglsLOW = lambda z : -0.182 * z + 41.123
psLOW = lambda z: 0.741 * z - 2.534
saundersFctLow=lambda logl, z : 10**psLOW(z) * (10**logl/10**loglsLOW(z))**(a+1) * n.e**( -n.log10( 1 +10**logl/10**loglsLOW(z))**2./(2*sig**2.))


#log_{10}(L_*) & =( -0.182 \pm  0.21 ) z + ( 41.272 \pm  0.149 ) \\
#log_{10}(\Phi_*) & =( 0.741 \pm  0.347 ) z + ( -2.78 \pm  0.246 ) \\
#\alpha &= -1.8 \pm 1.1 \\% ({\rm fixed}) \\
#\sigma &= 0.53 \; ({\rm fixed}) \\ 

logl, phi, phierr, z = [], [], [] ,[]
for ii in range(len(files)):
    x, y, ye = n.loadtxt(files[ii],unpack=True)
    logl.append(n.log10(x))
    phi.append(y)
    phierr.append(ye)
    zz = float(files[ii].split('_')[-1].split('.')[-2][1:])/1000.
    if zz>0.5 :
        z.append( n.ones_like(ye)*zz )
    if zz<0.5 :
        z.append( n.ones_like(ye)*1.2 )

logl = n.hstack((logl))
phi = n.hstack((phi))
phierr = n.hstack((phierr))
z = n.hstack(( z ))

ys = n.arange(z.min()-2*0.05,z.max()+2*0.05,0.05)
xs = n.arange(logl.min()-2*0.1, logl.max()+2*0.1,0.1)
X,Y = n.meshgrid(xs,ys)
Z = saundersFct(X,Y)
Z_up = saundersFctUp(X,Y)
Z_low = saundersFctLow(X,Y)


fig = p.figure(1,(9,9))
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(X, Y, n.log10(Z))#, rstride=10, cstride=10)
ax.plot_wireframe(X, Y, n.log10(Z_up))#, rstride=10, cstride=10)
ax.plot_wireframe(X, Y, n.log10(Z_low))#, rstride=10, cstride=10)

ok = (phi>phierr)

sc1 = ax.scatter(logl[ok],z[ok],n.log10(phi[ok]), s=n.ones_like(z[ok])*5, c='r', marker='o')
sc1.set_edgecolor('face')
#sc1 = ax.errorbar(logl[ok],z[ok],n.log10(phi[ok]), yerr=phierr[ok]/phi[ok])

ax.legend()
ax.set_xlabel(r'log ($L[O_{II}]/$ [erg s$^{-1}$])')
ax.set_ylabel('redshift')
ax.set_ylim((0.5,1.25))
ax.set_zlabel(r'$\Phi$ [Mpc$^{-3}/$dlog L]')
#ax.set_yscale('log')
#ax.set_zscale('log')
p.savefig(join(Pdir ,line, "LF-evolution.pdf"))
p.show()


fig = p.figure(1,(9,9))
ax = fig.add_subplot(111, projection='3d')

sc1 = ax.scatter(logl[ok],z[ok],phi[ok]/saundersFct(logl[ok],z[ok]), s=n.ones_like(z[ok])*3, c='r', marker='o', rasterized=True)
sc1.set_edgecolor('face')

ax.legend()
ax.set_xlabel(r'log ($L[O_{II}]/$ [erg s$^{-1}$])')
ax.set_ylabel('redshift')
ax.set_ylim((0.5,1.25))
ax.set_zlim((0,2))
ax.set_zlabel(r'data / model')
#ax.set_yscale('log')
#ax.set_zscale('log')
p.savefig(join(Pdir ,line, "LF-evolution-data-model-ratio.pdf"))
p.show()




sys.exit()

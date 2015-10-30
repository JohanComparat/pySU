import numpy as n
import matplotlib.pyplot as p
import glob
import sys

Pdir = "/Volumes/data/BigMD/2.5Gpc/propertiesAtVir/plots/"
dir = "/Volumes/data/BigMD/2.5Gpc/propertiesAtVir/"

Rbins = n.loadtxt(dir+'Rvir.bins')
Vbins = n.loadtxt(dir+'Vmax.bins')
Mbins = n.loadtxt(dir+'Mvir.bins')

# mvir rvir plane
centralList = n.array(glob.glob(dir+"hist2d-central-Mvir-Rvir-?.?????.dat"))
centralList.sort()
print centralList

X,Y= n.meshgrid(Rbins, Mbins)

volume = (2500.)**3. 
norm = volume * n.median(Mbins[1:]-Mbins[:-1]) *n.median(Rbins[1:]-Rbins[:-1]) 
print norm , n.log10(0.8/norm)

p.figure(1,(10,6))
for iii,el in enumerate(centralList) :
    z=1./float(el.split('-')[-1][:-4])-1.
    data = n.loadtxt(el,unpack=True)
    p.title("z="+str(n.round(z,3)))
    p.pcolormesh(X, Y, n.log10(data.T/norm),vmin=n.log10(0.8/norm),vmax=-4,rasterized=True)
    cb = p.colorbar(shrink=0.7)
    cb.set_label('N/(V dlogM dAR) ')
    p.ylim((11,16))
    p.ylabel(r'log $M_{vir}$ [km/s]')
    p.xlim((0,3.5))
    p.xlabel('Rvir [kpc/h]')
    p.savefig(Pdir + "plot-hist2d-central-Mvir-Rvir-"+str(n.round(z,5))+".pdf")
    p.clf()


# vmax rvir plane
centralList = n.array(glob.glob(dir+"hist2d-central-Vmax-Rvir-?.?????.dat"))
centralList.sort()
print centralList

X,Y= n.meshgrid(Rbins, Vbins)

volume = (2500.)**3. 
norm = volume * n.median(Vbins[1:]-Vbins[:-1]) *n.median(Rbins[1:]-Rbins[:-1]) 
p.figure(1,(10,6))
for iii,el in enumerate(centralList) :
    z=1./float(el.split('-')[-1][:-4])-1.
    data = n.loadtxt(el,unpack=True)
    p.title("z="+str(n.round(z,3)))
    p.pcolormesh(X, Y, n.log10(data.T/norm),vmin=n.log10(0.8/norm),vmax=-4,rasterized=True)
    cb = p.colorbar(shrink=0.7)
    cb.set_label('N/(V dlogV dAR) ')
    p.ylim((1.5,3.5))
    p.ylabel(r'log $V_{max}$ [km/s]')
    p.xlim((0,50000))
    p.xlabel('Mvir Msun/h/yr')
    p.savefig(Pdir + "plot-hist2d-central-Vmax-Rvir-"+str(n.round(z,5))+".pdf")
    p.clf()


# Mvir function
centralList = n.array(glob.glob(dir+"hist-sat-Mvir-?.?????.dat"))
centralList.sort()
print centralList

p.figure(1,(10,6))

for iii,el in enumerate(centralList) :
    z=1./float(el.split('-')[-1][:-4])-1.
    b0,b1,val = n.loadtxt(el,unpack=True)
    p.plot((10**b0+10**b1)/2.,val/(volume * (b1 - b0) ),label="z="+str(n.round(z,3)), rasterized=True)
    
    if iii%15==0 :
        p.axhline(1/(volume*n.median((b1 - b0))),label=r'1/(dlogMx2500$^3$)')
        p.axvline(23593750000.0*100, label='100 Mp')
        p.xlim((23593750000.0*50,5e16))
        p.ylim((0.5/(volume*n.median((b1 - b0))), 1e-1))
        p.xlabel(r'$M_{vir}$ [M$_\odot$/h]')
        p.ylabel(r'N($V_1<M_{vir}<V_2$) / dlog Mvir / Volume [ h/M$_\odot$ . h3/Mpc3 ]')
        p.xscale('log')    
        p.yscale('log')    
        p.legend()
        p.grid()
        p.savefig(Pdir + "plot-hist-sat-Mvir-"+str(iii)+".pdf")
        p.clf()


# Rvir function
centralList = n.array(glob.glob(dir+"hist-sat-Rvir-?.?????.dat"))
centralList.sort()
print centralList

p.figure(1,(10,6))

for iii,el in enumerate(centralList) :
    z=1./float(el.split('-')[-1][:-4])-1.
    b0,b1,val = n.loadtxt(el,unpack=True)
    p.plot((10**b0+10**b1)/2.,val/(volume * (b1 - b0) ),label="z="+str(n.round(z,3)), rasterized=True)
    
    if iii%15==0 :
        p.axhline(1/(volume*n.median((b1 - b0))),label=r'1/(dlogMx2500$^3$)')
        p.axvline(23593750000.0*100, label='100 Mp')
        p.xlim((23593750000.0*50,5e16))
        p.ylim((0.5/(volume*n.median((b1 - b0))), 1e-1))
        p.xlabel('Rvir [kpc/h]')
        p.ylabel(r'N($V_1<R_{vir}<V_2$) / dlog Rvir / Volume [ h/M$_\odot$ . h3/Mpc3 ]')
        p.xscale('log')    
        p.yscale('log')    
        p.legend()
        p.grid()
        p.savefig(Pdir + "plot-hist-sat-Rvir-"+str(iii)+".pdf")
        p.clf()


centralList = n.array(glob.glob(dir+"hist-central-Vmax-?.?????.dat"))
centralList.sort()

volume = (2500.)**3. 

p.figure(1,(10,6))

for iii,el in enumerate(centralList) :
    z=1./float(el.split('-')[-1][:-4])-1.
    b0,b1,val = n.loadtxt(el,unpack=True)
    p.plot((10**b0+10**b1)/2.,val/(volume * (b1 - b0) ),label="z="+str(n.round(z,3)), rasterized=True)

    if iii%15==0 :
        p.axhline(1/(volume*n.median((b1 - b0))),label=r'1/(0.01x2500$^3$)')
        p.xlim((100,1e5))
        p.ylim((0.5/(volume*n.median((b1 - b0))), 1e-1))
        p.xlabel(r'$V_{max}$ [km/s]')
        p.ylabel(r'N($V_1<V_{max}<V_2$) / dlog Vmax / Volume [ s/km . h3/Mpc3 ]')
        p.xscale('log')    
        p.yscale('log')    
        p.legend()
        p.grid()
        p.savefig(Pdir + "plot-hist-central-Vmax-"+str(iii)+".pdf")
        p.clf()


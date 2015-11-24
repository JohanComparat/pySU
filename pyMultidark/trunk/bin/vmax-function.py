from scipy.interpolate import interp1d
"""

ls Multidark-lightcones/MD_*/properties/vmax/hist-Central*1.0*.dat
ls Multidark-lightcones/MD_*/properties/vmax/hist-Central*0.9*.dat
ls Multidark-lightcones/MD_*/properties/vmax/hist-Central*0.8*.dat
ls Multidark-lightcones/MD_*/properties/vmax/hist-Central*0.7*.dat
ls Multidark-lightcones/MD_*/properties/vmax/hist-Central*0.6*.dat
ls Multidark-lightcones/MD_*/properties/vmax/hist-Central*0.5*.dat
ls Multidark-lightcones/MD_*/properties/vmax/hist-Central*0.4*.dat
ls Multidark-lightcones/MD_*/properties/vmax/hist-Central*0.3*.dat
ls Multidark-lightcones/MD_*/properties/vmax/hist-Central*0.2*.dat
"""
import numpy as n
import matplotlib
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p
import glob
import sys
from scipy.optimize import curve_fit
import cPickle

dir = "/Volumes/data/BigMD/vmaxFunction/"
Pdir = "/Volumes/data/BigMD/vmaxFunction/plots/"

def getVF(b0, b1, val,volume,label="SMD",completeness = 100, maxV=1500,errFactor=1.):
    Nc = n.array([ n.sum(val[ii:]) for ii in range(len(val)) ])
    xData_A = (10**b0+10**b1)/2.
    yData_A = Nc/(volume ) 
    yDataErr_A = errFactor/val**0.5
    sel = (yDataErr_A!=n.inf)&(yData_A>0) &(yData_A * volume >= 1) &(xData_A > completeness)&(xData_A < maxV)
    xData = xData_A[sel]
    yData = yData_A[sel]
    yDataErr = yDataErr_A[sel]*yData_A[sel]
    return xData,yData,yDataErr,volume


def plotVFv3(b0, b1, val,volume,label="SMD",errFactor=1.):
    Nc = n.array([ n.sum(val[ii:]) for ii in range(len(val)) ])
    xData_A = (10**b0+10**b1)/2.
    yData_A = Nc/(volume ) 
    yDataErr_A = errFactor/val**0.5
    sel = (yDataErr_A!=n.inf)&(yData_A>0) &(yData_A * volume >= 1)
    xData = xData_A[sel]
    yData = yData_A[sel]
    yDataErr = yDataErr_A[sel]*yData_A[sel]
    p.errorbar(xData,xData**3.1*yData,yerr=xData**3.1*yDataErr,elinewidth=1, label=label, rasterized=True)
    #print yData
    #p.axhline(1/(volume),label=r'1/V '+label,ls='dashed')

def plot_MVIR_velocity(matrix,xname,yname,label):



vf = lambda v, A, v0, alpha, beta : 10**A * v**beta * n.e**(- (v/10**v0)**alpha )
vx = n.logspace(1.7,3.4,100)
pl = vf(vx,4,3,2.15,-2.9)


######### z = 0.00 a = 1.0  hist2D Mvir Vmax ############

SMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_0.4Gpc/properties/vmax-mvir/hist2d-Central-vmax-mvir-1.00000.dat"
mvirBinsFile = "/Volumes/data/BigMD/Multidark-lightcones/MD_0.4Gpc/properties/vmax-mvir/mvir.bins"
vmaxBinsFile = "/Volumes/data/BigMD/Multidark-lightcones/MD_0.4Gpc/properties/vmax-mvir/vmax.bins"

MDPLfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_1Gpc_new_rockS/properties/vmax-mvir/hist2d-Central-vmax-mvir-1.00000.dat"
BigMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_2.5Gpc/properties/vmax-mvir/hist2d-Central-vmax-mvir-1.00000.dat"
HMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_4Gpc/properties/vmax-mvir/hist2d-Central-vmax-mvir-1.00000.dat"


H_04 = n.loadtxt(SMDfile)
mvir_bins_04 = n.loadtxt(mvirBinsFile)
vmax_bins_04 = n.loadtxt(vmaxBinsFile)
X,Y = n.meshgrid(( mvir_bins_04[1:] + mvir_bins_04[:-1])/2., (vmax_bins_04[1:]+vmax_bins_04[:-1])/2.)

sel = (H_04>0)

n.sum(H_04, axis=0)[n.sum(H_04, axis=0)>0]
n.sum(H_04, axis=1)[n.sum(H_04, axis=1)>0]

relation = interp1d(( mvir_bins_04[1:] + mvir_bins_04[:-1])/2.,H_04[180])
H_04[180].sum()



p.figure(1,(6,6))
p.axes([0.17,0.17,0.75,0.75])
H_04 = n.loadtxt(SMDfile)
p.contour(X,Y,n.log10(H_04))#, colors='r')
H_10 = n.loadtxt(MDPLfile)
p.contour(X,Y,n.log10(H_10))#, colors='r')
#H_25 = n.loadtxt(BigMDfile)
#p.contour(X,Y,H_25)
H_40 = n.loadtxt(HMDfile)
p.contour(X,Y,n.log10(H_40))#, colors='r')
p.xlabel(r'M_${vir}$')
p.ylabel(r'V_${max}$')
p.colorbar()
p.savefig(Pdir + "mvir-vmax-z0.00.pdf")
p.show()

######### z = 0.00 a = 1.0 ############

SMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_0.4Gpc/properties/vmax-mvir/hist-Central-vmax-1.00000.dat"
MDPLfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_1Gpc_new_rockS/properties/vmax-mvir/hist-Central-vmax-1.00000.dat"
BigMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_2.5Gpc/properties/vmax-mvir/hist-Central-vmax-1.00000.dat"
HMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_4Gpc/properties/vmax-mvir/hist-Central-vmax-1.00000.dat"


p.figure(1,(6,6))
p.axes([0.17,0.17,0.75,0.75])

b0_04, b1_04, val_04 = n.loadtxt(SMDfile,unpack=True)
z_04 = 1/float(SMDfile.split('-')[-1][:-4])-1
#xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,label="SMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,completeness = 50, maxV = 150, errFactor=10.)
p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=2,fmt='none', label="SMD z="+str(n.round(z_04,3)), rasterized=True)

b0_10, b1_10, val_10 = n.loadtxt(MDPLfile,unpack=True)
z_10 = 1/float(MDPLfile.split('-')[-1][:-4])-1
#xData_10,yData_10,yDataErr_10,volume_10 = getVF(b0_10, b1_10, val_10, 1000.**3.,label="MDPL",completeness = 0, maxV = 5000)
#p.errorbar(xData_10,yData_10,yerr=yDataErr_10,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_10,yData_10,yDataErr_10,volume_10 = getVF(b0_10, b1_10, val_10, 1000.**3.,completeness = 151, maxV = 300, errFactor=10.)
p.errorbar(xData_10,yData_10,yerr=yDataErr_10,elinewidth=2,fmt='none', label="MDPL z="+str(n.round(z_10,3)), rasterized=True)

b0_25, b1_25, val_25 = n.loadtxt(BigMDfile,unpack=True)
z_25 = 1/float(BigMDfile.split('-')[-1][:-4])-1
#xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,label="BigMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,completeness = 301, maxV = 500)
p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=2,fmt='none', label="BigMD z="+str(n.round(z_25,3)), rasterized=True)

b0_40, b1_40, val_40 = n.loadtxt(HMDfile,unpack=True)
z_40 = 1/float(HMDfile.split('-')[-1][:-4])-1
#xData_40,yData_40,yDataErr_40,volume_40 = getVF(b0_40, b1_40, val_40,4000.**3.,label="HMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_40,yData_40,yerr=yDataErr_40,elinewidth=1, ecolor='k',fmt='none', rasterized=True)
xData_40,yData_40,yDataErr_40,volume_40 = getVF(b0_40, b1_40, val_40,4000.**3.,completeness = 501, maxV = 2500)
p.errorbar(xData_40,yData_40,yerr=yDataErr_40,elinewidth=2,fmt='none', label="HMD z="+str(n.round(z_40,3)), rasterized=True)
p.axhline(1/(volume_40),label=r'1/V(HMD)',ls='dashed')

xData = n.hstack((xData_04, xData_25, xData_40))
yData = n.hstack((yData_04, yData_25, yData_40))
yDataErr = n.hstack((yDataErr_04, yDataErr_25, yDataErr_40))
res, cov = curve_fit(vf, xData, yData, sigma = yDataErr, p0 = (4,3,2.15,-3.2) , maxfev = 5000000)

p.plot(vx, vf(vx,res[0], res[1], res[2], res[3]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))

p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$)  [ h$^3$ Mpc$^{-3}$ ]')
p.xscale('log')
p.xlim((40,3000))
p.ylim((0.5/(volume_40),1))
p.yscale('log')
p.legend(loc=3)
p.grid()
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.savefig(Pdir + "fit-perBox-z0.00.pdf")
p.clf()


p.figure(2,(6,6))
p.axes([0.17,0.15,0.75,0.75])
p.fill_between(xData,1+yDataErr / yData,1-yDataErr / yData,rasterized=True,alpha=0.5)
p.plot(xData, yData / vf(xData,res[0], res[1], res[2], res[3]),'k')
p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$) / best fit')
p.xscale('log')
p.ylim((0.9,1.1))
p.xlim((40,3000))
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.grid()
p.savefig(Pdir + "fit-ratio-z0.00.pdf")
p.clf()

print n.round(z_25,3), " & ", n.round(res,3), "\\\\"



######### z = 0.05 a = 0.950 ############

SMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_0.4Gpc/properties/vmax/hist-Central-vmax-0.95600.dat"
MDPLfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_1Gpc_new_rockS/properties/vmax/hist-Central-vmax-0.95670.dat"
BigMDfile = "/Volumes/data/BigMD/2.5Gpc/vmax/hist-Central-vmax-0.95600.dat"
#HMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_4Gpc/properties/vmax/hist-Central-vmax-0.25320.dat"


p.figure(1,(6,6))
p.axes([0.17,0.17,0.75,0.75])

b0_25, b1_25, val_25 = n.loadtxt(BigMDfile,unpack=True)
z_25 = 1/float(BigMDfile.split('-')[-1][:-4])-1
#xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,label="BigMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,completeness = 220, maxV = 2500)
p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=2,fmt='none', label="BigMD z="+str(n.round(z_25,3)), rasterized=True)
p.axhline(1/(volume_25),label=r'1/V(BigMD)',ls='dashed')

b0_10, b1_10, val_10 = n.loadtxt(MDPLfile,unpack=True)
z_10 = 1/float(MDPLfile.split('-')[-1][:-4])-1
#xData_10,yData_10,yDataErr_10,volume_10 = getVF(b0_10, b1_10, val_10, 1000.**3.,label="MDPL",completeness = 0, maxV = 5000)
#p.errorbar(xData_10,yData_10,yerr=yDataErr_10,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_10,yData_10,yDataErr_10,volume_10 = getVF(b0_10, b1_10, val_10, 1000.**3.,completeness = 120, maxV = 220)
p.errorbar(xData_10,yData_10,yerr=yDataErr_10,elinewidth=2,fmt='none', label="MDPL z="+str(n.round(z_10,3)), rasterized=True)

b0_04, b1_04, val_04 = n.loadtxt(SMDfile,unpack=True)
z_04 = 1/float(SMDfile.split('-')[-1][:-4])-1
#xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,label="SMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,completeness = 50, maxV = 120)
p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=2,fmt='none', label="SMD z="+str(n.round(z_04,3)), rasterized=True)

xData = n.hstack((xData_04, xData_10, xData_25))
yData = n.hstack((yData_04, yData_10, yData_25))
yDataErr = n.hstack((yDataErr_04, yDataErr_10, yDataErr_25))
res, cov = curve_fit(vf, xData, yData, sigma = yDataErr, p0 = (4,3,2.15,-3.2) , maxfev = 5000000)

p.plot(vx, vf(vx,res[0], res[1], res[2], res[3]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))

p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$)  [ h$^3$ Mpc$^{-3}$ ]')
p.xscale('log')
p.yscale('log')
p.ylim((0.5/(volume_40),1))
p.xlim((40,3000))
p.legend(loc=3)
p.grid()
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.savefig(Pdir + "fit-perBox-z0.05.pdf")
p.clf()


p.figure(2,(6,6))
p.axes([0.17,0.15,0.75,0.75])
p.fill_between(xData,1+yDataErr / yData,1-yDataErr / yData,rasterized=True,alpha=0.5)
p.plot(xData, yData / vf(xData,res[0], res[1], res[2], res[3]),'k')
p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$) / best fit')
p.xscale('log')
p.ylim((0.95,1.05))
p.xlim((40,3000))
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.grid()
p.savefig(Pdir + "fit-ratio-z0.05.pdf")
p.clf()

print n.round(z_25,3), " & ", n.round(res,3), "\\\\"



######### z = 0.42 a = 0.70 ############

SMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_0.4Gpc/properties/vmax/hist-Central-vmax-0.70030.dat"
MDPLfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_1Gpc_new_rockS/properties/vmax/hist-Central-vmax-0.70160.dat"
BigMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_2.5Gpc/properties/vmax/hist-Central-vmax-0.70030.dat"
#HMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_4Gpc/properties/vmax/hist-Central-vmax-0.25320.dat"


p.figure(1,(6,6))
p.axes([0.17,0.17,0.75,0.75])

b0_25, b1_25, val_25 = n.loadtxt(BigMDfile,unpack=True)
z_25 = 1/float(BigMDfile.split('-')[-1][:-4])-1
#xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,label="BigMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,completeness = 220, maxV = 2500)
p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=2,fmt='none', label="BigMD z="+str(n.round(z_25,3)), rasterized=True)
p.axhline(1/(volume_25),label=r'1/V(BigMD)',ls='dashed')

b0_10, b1_10, val_10 = n.loadtxt(MDPLfile,unpack=True)
z_10 = 1/float(MDPLfile.split('-')[-1][:-4])-1
#xData_10,yData_10,yDataErr_10,volume_10 = getVF(b0_10, b1_10, val_10, 1000.**3.,label="MDPL",completeness = 0, maxV = 5000)
#p.errorbar(xData_10,yData_10,yerr=yDataErr_10,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_10,yData_10,yDataErr_10,volume_10 = getVF(b0_10, b1_10, val_10, 1000.**3.,completeness = 120, maxV = 220)
p.errorbar(xData_10,yData_10,yerr=yDataErr_10,elinewidth=2,fmt='none', label="MDPL z="+str(n.round(z_10,3)), rasterized=True)

b0_04, b1_04, val_04 = n.loadtxt(SMDfile,unpack=True)
z_04 = 1/float(SMDfile.split('-')[-1][:-4])-1
#xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,label="SMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,completeness = 50, maxV = 120)
p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=2,fmt='none', label="SMD z="+str(n.round(z_04,3)), rasterized=True)

xData = n.hstack((xData_04, xData_10, xData_25))
yData = n.hstack((yData_04, yData_10, yData_25))
yDataErr = n.hstack((yDataErr_04, yDataErr_10, yDataErr_25))
res, cov = curve_fit(vf, xData, yData, sigma = yDataErr, p0 = (4,3,2.15,-3.2) , maxfev = 5000000)

p.plot(vx, vf(vx,res[0], res[1], res[2], res[3]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))

p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$)  [ h$^3$ Mpc$^{-3}$ ]')
p.xscale('log')
p.yscale('log')
p.ylim((0.5/(volume_40),1))
p.xlim((40,3000))
p.legend(loc=3)
p.grid()
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.savefig(Pdir + "fit-perBox-z0.42.pdf")
p.clf()


p.figure(2,(6,6))
p.axes([0.17,0.15,0.75,0.75])
p.fill_between(xData,1+yDataErr / yData,1-yDataErr / yData,rasterized=True,alpha=0.5)
p.plot(xData, yData / vf(xData,res[0], res[1], res[2], res[3]),'k')
p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$) / best fit')
p.xscale('log')
p.ylim((0.95,1.05))
p.xlim((40,3000))
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.grid()
p.savefig(Pdir + "fit-ratio-z0.42.pdf")
p.clf()

print n.round(z_25,3), " & ", n.round(res,3), "\\\\"


######### z = 0.88 a = 0.5 ############

SMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_0.4Gpc/properties/vmax/hist-Central-vmax-0.53000.dat"
BigMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_2.5Gpc/properties/vmax/hist-Central-vmax-0.53000.dat"
HMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_4Gpc/properties/vmax/hist-Central-vmax-0.53780.dat"


p.figure(1,(6,6))
p.axes([0.17,0.17,0.75,0.75])
b0_40, b1_40, val_40 = n.loadtxt(HMDfile,unpack=True)
z_40 = 1/float(HMDfile.split('-')[-1][:-4])-1
#xData_40,yData_40,yDataErr_40,volume_40 = getVF(b0_40, b1_40, val_40,4000.**3.,label="HMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_40,yData_40,yerr=yDataErr_40,elinewidth=1, ecolor='k',fmt='none', rasterized=True)
xData_40,yData_40,yDataErr_40,volume_40 = getVF(b0_40, b1_40, val_40,4000.**3.,completeness = 1000, maxV = 2500)
p.errorbar(xData_40,yData_40,yerr=yDataErr_40,elinewidth=2,fmt='none', label="HMD z="+str(n.round(z_40,3)), rasterized=True)
p.axhline(1/(volume_40),label=r'1/V(HMD)',ls='dashed')

b0_25, b1_25, val_25 = n.loadtxt(BigMDfile,unpack=True)
z_25 = 1/float(BigMDfile.split('-')[-1][:-4])-1
#xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,label="BigMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,completeness = 200, maxV = 1000)
p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=2,fmt='none', label="BigMD z="+str(n.round(z_25,3)), rasterized=True)

b0_04, b1_04, val_04 = n.loadtxt(SMDfile,unpack=True)
z_04 = 1/float(SMDfile.split('-')[-1][:-4])-1
#xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,label="SMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,completeness = 50, maxV = 200)
p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=2,fmt='none', label="SMD z="+str(n.round(z_04,3)), rasterized=True)

xData = n.hstack((xData_04, xData_25, xData_40))
yData = n.hstack((yData_04, yData_25, yData_40))
yDataErr = n.hstack((yDataErr_04, yDataErr_25, yDataErr_40))
res, cov = curve_fit(vf, xData, yData, sigma = yDataErr, p0 = (4,3,2.15,-3.2) , maxfev = 5000000)

p.plot(vx, vf(vx,res[0], res[1], res[2], res[3]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))

p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$)  [ h$^3$ Mpc$^{-3}$ ]')
p.xscale('log')
p.xlim((40,3000))
p.ylim((0.5/(volume_40),1))
p.yscale('log')
p.legend(loc=3)
p.grid()
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.savefig(Pdir + "fit-perBox-z0.88.pdf")
p.clf()


p.figure(2,(6,6))
p.axes([0.17,0.15,0.75,0.75])
p.fill_between(xData,1+yDataErr / yData,1-yDataErr / yData,rasterized=True,alpha=0.5)
p.plot(xData, yData / vf(xData,res[0], res[1], res[2], res[3]),'k')
p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$) / best fit')
p.xscale('log')
p.ylim((0.95,1.05))
p.xlim((40,3000))
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.grid()
p.savefig(Pdir + "fit-ratio-z0.88.pdf")
p.clf()


print n.round(z_25,3), " & ", n.round(res,3), "\\\\"



######### z = 3 a = 0.25 ############

SMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_0.4Gpc/properties/vmax/hist-Central-vmax-0.24800.dat"
MDPLfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_1Gpc_new_rockS/properties/vmax/hist-Central-vmax-0.24770.dat"
BigMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_2.5Gpc/properties/vmax/hist-Central-vmax-0.25700.dat"
HMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_4Gpc/properties/vmax/hist-Central-vmax-0.25320.dat"


p.figure(1,(6,6))
p.axes([0.17,0.17,0.75,0.75])
b0_40, b1_40, val_40 = n.loadtxt(HMDfile,unpack=True)
z_40 = 1/float(HMDfile.split('-')[-1][:-4])-1
#xData_40,yData_40,yDataErr_40,volume_40 = getVF(b0_40, b1_40, val_40,4000.**3.,label="HMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_40,yData_40,yerr=yDataErr_40,elinewidth=1, ecolor='k',fmt='none', rasterized=True)
xData_40,yData_40,yDataErr_40,volume_40 = getVF(b0_40, b1_40, val_40,4000.**3.,completeness = 1030, maxV = 2500)
p.errorbar(xData_40,yData_40,yerr=yDataErr_40,elinewidth=2,fmt='none', label="HMD z="+str(n.round(z_40,3)), rasterized=True)
p.axhline(1/(volume_40),label=r'1/V(HMD)',ls='dashed')

b0_25, b1_25, val_25 = n.loadtxt(BigMDfile,unpack=True)
z_25 = 1/float(BigMDfile.split('-')[-1][:-4])-1
#xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,label="BigMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,completeness = 220, maxV = 1030)
p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=2,fmt='none', label="BigMD z="+str(n.round(z_25,3)), rasterized=True)

b0_10, b1_10, val_10 = n.loadtxt(MDPLfile,unpack=True)
z_10 = 1/float(MDPLfile.split('-')[-1][:-4])-1
#xData_10,yData_10,yDataErr_10,volume_10 = getVF(b0_10, b1_10, val_10, 1000.**3.,label="MDPL",completeness = 0, maxV = 5000)
#p.errorbar(xData_10,yData_10,yerr=yDataErr_10,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_10,yData_10,yDataErr_10,volume_10 = getVF(b0_10, b1_10, val_10, 1000.**3.,completeness = 120, maxV = 220)
p.errorbar(xData_10,yData_10,yerr=yDataErr_10,elinewidth=2,fmt='none', label="MDPL z="+str(n.round(z_10,3)), rasterized=True)

b0_04, b1_04, val_04 = n.loadtxt(SMDfile,unpack=True)
z_04 = 1/float(SMDfile.split('-')[-1][:-4])-1
#xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,label="SMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,completeness = 50, maxV = 120)
p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=2,fmt='none', label="SMD z="+str(n.round(z_04,3)), rasterized=True)

xData = n.hstack((xData_04, xData_10, xData_25, xData_40))
yData = n.hstack((yData_04, yData_10, yData_25, yData_40))
yDataErr = n.hstack((yDataErr_04, yDataErr_10, yDataErr_25, yDataErr_40))
res, cov = curve_fit(vf, xData, yData, sigma = yDataErr, p0 = (4,3,2.15,-3.2) , maxfev = 5000000)

p.plot(vx, vf(vx,res[0], res[1], res[2], res[3]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))

p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$)  [ h$^3$ Mpc$^{-3}$ ]')
p.xscale('log')
p.yscale('log')
p.ylim((0.5/(volume_40),1))
p.xlim((40,3000))
p.legend(loc=3)
p.grid()
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.savefig(Pdir + "fit-perBox-z3.0.pdf")
p.clf()


p.figure(2,(6,6))
p.axes([0.17,0.15,0.75,0.75])
p.fill_between(xData,1+yDataErr / yData,1-yDataErr / yData,rasterized=True,alpha=0.5)
p.plot(xData, yData / vf(xData,res[0], res[1], res[2], res[3]),'k')
p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$) / best fit')
p.xscale('log')
p.ylim((0.9,1.1))
p.xlim((40,3000))
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.grid()
p.savefig(Pdir + "fit-ratio-z3.0.pdf")
p.clf()


print n.round(z_25,3), " & ", n.round(res,3), "\\\\"





sys.exit()















p.figure(1,(6,6))
p.axes([0.2,0.2,0.75,0.75])

b0_10, b1_10, val_10 = n.loadtxt("Multidark-lightcones/MD_1Gpc_new_rockS/properties/vmax/hist-Central-vmax-0.24770.dat",unpack=True)
xData,yData,yDataErr,volume = plotVF(b0_10, b1_10, val_10,1000.**3.,label="MDPL",completeness = 100, maxV = 1500)

p.errorbar(xData,yData,yerr=yDataErr,elinewidth=1, label=label, rasterized=True)
p.axhline(1/(volume),label=r'1/V '+label,ls='dashed')

b0_25, b1_25, val_25 = n.loadtxt("Multidark-lightcones/MD_2.5Gpc/properties/vmax/hist-Central-vmax-0.25700.dat",unpack=True)
plotVF(b0_25, b1_25, val_25, 2500.**3.,label="BigMD")


plotVF(b0_25, b1_25, val_25, 4000.**3.,label="HMD")

p.xlim((50,4000))
p.ylim((0.5/(4000.**3.), 1))
p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$)  [ h$^3$ Mpc$^{-3}$ ]')
p.xscale('log')
p.yscale('log')
p.legend(loc=3)
p.savefig(Pdir + "VmaxF-cumulative-central-z3.pdf")


p.figure(2,(6,6))
p.axes([0.2,0.2,0.75,0.75])

b0_04, b1_04, val_04 = n.loadtxt("Multidark-lightcones/MD_0.4Gpc/properties/vmax/hist-Central-vmax-0.71830.dat",unpack=True)
plotVF(b0_04, b1_04, val_04,400.**3.,label="SMD")

b0_10, b1_10, val_10 = n.loadtxt("Multidark-lightcones/MD_1Gpc_new_rockS/properties/vmax/hist-Central-vmax-0.71730.dat",unpack=True)
plotVF(b0_10, b1_10, val_10,1000.**3.,label="MDPL")

b0_25, b1_25, val_25 = n.loadtxt("2.5Gpc/vmax/hist-Central-vmax-0.71830.dat",unpack=True)
plotVF(b0_25, b1_25, val_25, 2500.**3.,label="BigMD")


p.xlim((5,4000))
p.ylim((0.5/(2500.**3.), 10))
p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$)  [ h$^3$ Mpc$^{-3}$ ]')
p.xscale('log')
p.yscale('log')
p.legend(loc=3)
p.grid()
p.savefig(Pdir + "VmaxF-cumulative-central-z0.42.pdf")
p.clf()



p.figure(2,(6,6))
p.axes([0.2,0.2,0.75,0.75])

#b0_04, b1_04, val_04 = n.loadtxt("Multidark-lightcones/MD_0.4Gpc/properties/vmax/hist-Central-vmax-0.71830.dat",unpack=True)
#plotVFv3(b0_04, b1_04, val_04,400.**3.,label="SMD")

b0_10, b1_10, val_10 = n.loadtxt("Multidark-lightcones/MD_1Gpc_new_rockS/properties/vmax/hist-Central-vmax-0.71730.dat",unpack=True)
plotVFv3(b0_10, b1_10, val_10,1000.**3.,label="MDPL")

b0_25, b1_25, val_25 = n.loadtxt("2.5Gpc/vmax/hist-Central-vmax-0.71830.dat",unpack=True)
plotVFv3(b0_25, b1_25, val_25, 2500.**3.,label="BigMD")


p.xlim((100,400))
p.ylim((40000,50000))
p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'$V_{max}^{3.1}$N($>V_{max}$)  [ h$^3$ Mpc$^{-3}$ ]')
#p.xscale('log')
#p.yscale('log')
p.legend(loc=3)
p.grid()
p.savefig(Pdir + "VmaxFv3-cumulative-central-z0.42-zoom-o2.pdf")
p.clf()





centralList = n.array(glob.glob(dir+"hist-central-Vpeak-?.?????.dat"))
centralList.sort()

volume = (2500.)**3. 

# evolution plot

ids = [2,3,4,5,7,9,11,20,56]
p.figure(1,(6,6))
p.axes([0.2,0.2,0.75,0.75])
for iii in ids :
    el = centralList[iii]
    z=1./float(el.split('-')[-1][:-4])-1.
    b0,b1,val = n.loadtxt(el,unpack=True)
    Nc = n.array([ n.sum(val[ii:]) for ii in range(len(val)) ])
    xData_A = (10**b0+10**b1)/2.
    yData_A = Nc/(volume ) 
    yDataErr_A = 1/val**0.5
    sel = (yDataErr_A!=n.inf)&(yData_A>0) & (xData_A>20)&(yData_A * volume >= 1)
    xData = xData_A[sel]
    yData = yData_A[sel]
    yDataErr = yDataErr_A[sel]*yData_A[sel]
    p.errorbar(xData,yData,yerr=yDataErr,elinewidth=1, label="z="+str(n.round(z,3)), rasterized=True)

p.axhline(1/(volume),label=r'1/V',color='k',ls='dashed')
p.xlim((50,4000))
p.ylim((0.5/(volume), 1e-2))
p.xlabel(r'$V_{peak}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{peak}$)  [ h$^3$ Mpc$^{-3}$ ]')
p.xscale('log')
p.yscale('log')
p.legend(loc=3)
p.grid()
p.savefig(Pdir + "VpeakF-cumulative-central-z8-z0-evolution.pdf")
p.clf()

sys.exit()


filename = n.array(glob.glob("0.4Gpc/*/hist*Central-Vpeak-0.22*.dat"))[0]
def getData(filename,Vcut=0,volume = 400.**3.):
    zz=1./float(filename.split('-')[-1][:-4])-1.
    b0,b1,val = n.loadtxt(filename,unpack=True)
    Nc = n.array([ n.sum(val[ii:]) for ii in range(len(val)) ])
    xData_A = (10**b0+10**b1)/2.
    yData_A = Nc/(volume ) 
    yDataErr_A = 1/val**0.5
    sel = (yDataErr_A!=n.inf)&(yData_A>0) & (xData_A>Vcut)&(yData_A * volume >= 1)
    velocity = xData_A[sel]
    redshift = n.ones_like(velocity)*zz
    yData = yData_A[sel]
    yDataErr = yDataErr_A[sel]*yData_A[sel]
    return redshift, velocity, yData,yDataErr

z,v,y,yerr = getData(hist04)


hist10 = n.array(glob.glob("1Gpc/*/hist*Central-Vpeak-0.22*.dat"))[0]
z10=1./float(hist10.split('-')[-1][:-4])-1.
hist25 = n.array(glob.glob("2.5Gpc/*/hist*central-Vpeak-0.22*.dat"))[0]
z25=1./float(hist25.split('-')[-1][:-4])-1.



# cumulative velocity function

results, errors, redshifts, chi2Rs = [], [], [], []

for iii,el in enumerate(centralList) :
    z=1./float(el.split('-')[-1][:-4])-1.
    b0,b1,val = n.loadtxt(el,unpack=True)
    Nc = n.array([ n.sum(val[ii:]) for ii in range(len(val)) ])
    xData_A = (10**b0+10**b1)/2.
    yData_A = Nc/(volume ) 
    yDataErr_A = 1/val**0.5
    sel = (yDataErr_A!=n.inf)&(yData_A>0) & (xData_A>200)&(yData_A * volume >= 1)
    xData = xData_A[sel]
    yData = yData_A[sel]
    yDataErr = yDataErr_A[sel]*yData_A[sel]
    res, cov = curve_fit(vf, xData, yData, sigma = yDataErr, p0 = (4,3,2.15,-3.2) , maxfev = 5000000)
    chi2red = n.round(n.sum( (vf(xData,res[0], res[1], res[2], res[3])- yData)**2. / (yDataErr**2)) / (len(xData) - (4-1)),2)
    chi2Rs.append(chi2red)
    results.append(res)
    errors.append(cov)
    redshifts.append(z)
    p.figure(1,(6,6))
    p.axes([0.2,0.2,0.75,0.75])
    p.plot(vx, vf(vx,res[0], res[1], res[2], res[3]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))
    p.errorbar(xData,yData,yerr=yDataErr,color='k',elinewidth=2, label="z="+str(n.round(z,3)), rasterized=True)
    p.axhline(1/(volume),label=r'1/V',color='k',ls='dashed')
    p.xlim((100,4000))
    p.ylim((0.5/(volume), 1e-1))
    p.xlabel(r'$V_{peak}$ [km s$^{-1}$]')
    p.ylabel(r'N($>V_{peak}$)  [ h$^3$ Mpc$^{-3}$ ]')
    p.xscale('log')    
    p.yscale('log')    
    p.title(str(n.round(res[0],2))+" "+ str(n.round(res[1],2))+" "+ str(n.round(res[2],2))+" "+ str(n.round(res[3],2)) )
    p.legend()
    p.grid()
    p.savefig(Pdir + "VpeakF-cumulative-central-z-"+str(n.round(z,4))+".pdf")
    p.clf()

results = n.transpose(results)
errors = n.array(errors)
diag_err = n.array([ [el[0][0]**0.5, el[1][1]**0.5, el[2][2]**0.5, el[3][3]**0.5] for el in errors]).T
redshifts = n.array(redshifts)
chi2Rs = n.array(chi2Rs)

f=open(Pdir + "VpeakF-cumulative-central-z-param-fit-results.pkl",'w')
cPickle.dump([results, errors, diag_err, redshifts, chi2Rs],f)
f.close()

ok = (redshifts < 1.1)

p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshifts[ok],results[0][ok], diag_err[0][ok])
p.ylabel('log Amplitude')
p.xlabel('z')
p.grid()
p.savefig(Pdir + "VpeakF-cumulative-central-z-param-amplitude-evolution.pdf")
p.clf()
p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshifts[ok],results[1][ok], diag_err[1][ok])
p.ylabel(r'$V_0$')
p.xlabel('z')
p.grid()
p.savefig(Pdir + "VpeakF-cumulative-central-z-param-v0-evolution.pdf")
p.clf()
p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshifts[ok],results[2][ok], diag_err[2][ok])
p.ylabel(r'$\alpha$')
p.xlabel('z')
p.grid()
p.savefig(Pdir + "VpeakF-cumulative-central-z-param-alpha-evolution.pdf")
p.clf()
p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshifts[ok],results[3][ok], diag_err[3][ok])
p.ylabel(r'$\beta$')
p.xlabel('z')
p.grid()
p.savefig(Pdir + "VpeakF-cumulative-central-z-param-beta-evolution.pdf")
p.clf()
sys.exit()

Pdir = "/Volumes/data/BigMD/2.5Gpc/plots/"
dir = "/Volumes/data/BigMD/2.5Gpc/"

ARbins = n.loadtxt('/Volumes/data/BigMD/2.5Gpc/AccRate.bins')
Vbins = n.loadtxt('/Volumes/data/BigMD/2.5Gpc/Vpeak.bins')
Mbins = n.loadtxt('/Volumes/data/BigMD/2.5Gpc/Mpeak.bins')


centralList = n.array(glob.glob(dir+"hist2d-central-Mpeak-Accrate-?.?????.dat"))
centralList.sort()
print centralList

X,Y= n.meshgrid(ARbins, Mbins)

volume = (2500.)**3. 
norm = volume * n.median(Mbins[1:]-Mbins[:-1]) *n.median(ARbins[1:]-ARbins[:-1]) 
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
    p.ylabel(r'log $M_{peak}$ [km/s]')
    p.xlim((0,50000))
    p.xlabel('Accretion rate Msun/h/yr at Mpeak')
    p.savefig(Pdir + "plot-hist2d-central-Mpeak-Accrate-"+str(n.round(z,5))+".pdf")
    p.clf()


sys.exit()



centralList = n.array(glob.glob(dir+"hist2d-central-Vpeak-Accrate-?.?????.dat"))
centralList.sort()
print centralList

X,Y= n.meshgrid(ARbins, Vbins)

volume = (2500.)**3. 
norm = volume * n.median(Vbins[1:]-Vbins[:-1]) *n.median(ARbins[1:]-ARbins[:-1]) 
p.figure(1,(10,6))
for iii,el in enumerate(centralList) :
    z=1./float(el.split('-')[-1][:-4])-1.
    data = n.loadtxt(el,unpack=True)
    p.title("z="+str(n.round(z,3)))
    p.pcolormesh(X, Y, n.log10(data.T/norm),vmin=n.log10(0.8/norm),vmax=-4,rasterized=True)
    cb = p.colorbar(shrink=0.7)
    cb.set_label('N/(V dlogV dAR) ')
    p.ylim((1.5,3.5))
    p.ylabel(r'log $V_{peak}$ [km/s]')
    p.xlim((0,50000))
    p.xlabel('Growth Rate of Mpeak [z, z+0.5] Msun/h/yr')
    p.savefig(Pdir + "plot-hist2d-central-Vpeak-Accrate-"+str(n.round(z,5))+".pdf")
    p.clf()





centralList = n.array(glob.glob(dir+"hist-sat-Mpeak-?.?????.dat"))
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
        p.xlabel(r'$M_{peak}$ [M$_\odot$/h]')
        p.ylabel(r'N($V_1<M_{peak}<V_2$) / dlog Mpeak / Volume [ h/M$_\odot$ . h3/Mpc3 ]')
        p.xscale('log')    
        p.yscale('log')    
        p.legend()
        p.grid()
        p.savefig(Pdir + "plot-hist-sat-Mpeak-"+str(iii)+".pdf")
        p.clf()



centralList = n.array(glob.glob(dir+"hist-central-Mpeak-?.?????.dat"))
centralList.sort()
print centralList

volume = (2500.)**3. 

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
        p.xlabel(r'$M_{peak}$ [M$_\odot$/h]')
        p.ylabel(r'N($V_1<M_{peak}<V_2$) / dlog Mpeak / Volume [ h/M$_\odot$. h3/Mpc3 ]')
        p.xscale('log')    
        p.yscale('log')    
        p.legend()
        p.grid()
        p.savefig(Pdir + "plot-hist-central-Mpeak-"+str(iii)+".pdf")
        p.clf()


sys.exit()





volume = (2500.)**3. 

centralList = n.array(glob.glob(dir+"hist-sat-Vpeak-?.?????.dat"))
centralList.sort()


p.figure(1,(10,6))

for iii,el in enumerate(centralList) :
    z=1./float(el.split('-')[-1][:-4])-1.
    b0,b1,val = n.loadtxt(el,unpack=True)
    p.plot((10**b0+10**b1)/2.,val/(volume * (b1 - b0) ),label="z="+str(n.round(z,3)), rasterized=True)
    
    if iii%15==0 :
        p.axhline(1/(volume*n.median((b1 - b0))),label=r'1/(0.01x2500$^3$)')
        p.xlim((100,1e5))
        p.ylim((0.5/(volume*n.median((b1 - b0))), 1e-1))
        p.xlabel(r'$V_{peak}$ [km/s]')
        p.ylabel(r'N($V_1<V_{peak}<V_2$) / dlog Vpeak / Volume [ s/km . h3/Mpc3 ]')
        p.xscale('log')    
        p.yscale('log')    
        p.legend()
        p.grid()
        p.savefig(Pdir + "plot-hist-sat-Vpeak-"+str(iii)+".pdf")
        p.clf()



centralList = n.array(glob.glob(dir+"hist-central-Vpeak-?.?????.dat"))
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
        p.xlabel(r'$V_{peak}$ [km/s]')
        p.ylabel(r'N($V_1<V_{peak}<V_2$) / dlog Vpeak / Volume [ s/km . h3/Mpc3 ]')
        p.xscale('log')    
        p.yscale('log')    
        p.legend()
        p.grid()
        p.savefig(Pdir + "plot-hist-central-Vpeak-"+str(iii)+".pdf")
        p.clf()


sys.exit()




zzz=[]
for ii in range(len(centralList)):
    zzz.append( centralList[ii].split('-')[-1][:-4] )

Zset = list(set(zzz))

for el in Zset :
    centralList = n.array(glob.glob(dir+"hist-central-Vpeak-"+el+".dat"))
    centralList.sort()


    massB[:-1], massB[1:],  nnM.sum(axis=0)


    nnM = n.empty( [len(centralList),len(massB)-1] ) 
    nnV = n.empty( [len(centralList),len(vcirB)-1] )
    dataVC = n.empty( [len(centralList),len(vcirB)-1,len(concB)-1] )
    dataMC = n.empty( [len(centralList),len(massB)-1,len(concB)-1] )

    for jj in range(len(centralList)):
        f=open(centralList[jj],'r')
        nnMinter,nnVinter,nnCinter,dataMCinter,dataVCinter = cPickle.load(f)
        nnM[jj] = nnMinter
        nnV[jj] = nnVinter 
        dataMC[jj] = dataMCinter[0]
        dataVC[jj] = dataVCinter[0]
        f.close()


    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist-central-Mpeak-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",n.transpose([massB[:-1], massB[1:],  nnM.sum(axis=0)]))

    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist-central-Vpeak-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",n.transpose([vcirB[:-1], vcirB[1:],  nnV.sum(axis=0)]) )

    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist2d-central-Mpeak-Accrate-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",dataMC.sum(axis=0))

    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist2d-central-Vpeak-Accrate-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",dataVC.sum(axis=0))


    satList = n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/*"+box.get_snl()[ii].split('_')[-1][:-5]+"*MVAmatrixSatellite*.pkl"))
    satList.sort()

    nnM = n.empty( [len(satList),len(massB)-1] ) 
    nnV = n.empty( [len(satList),len(vcirB)-1] )
    dataVC = n.empty( [len(satList),len(vcirB)-1,len(concB)-1] )
    dataMC = n.empty( [len(satList),len(massB)-1,len(concB)-1] )

    for jj in range(len(satList)):
        f=open(satList[jj],'r')
        nnMinter,nnVinter,nnCinter,dataMCinter,dataVCinter = cPickle.load(f)
        nnM[jj] = nnMinter
        nnV[jj] = nnVinter 
        dataMC[jj] = dataMCinter[0]
        dataVC[jj] = dataVCinter[0]
        f.close()


    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist-sat-Mpeak-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",n.transpose([massB[:-1], massB[1:],  nnM.sum(axis=0)]))

    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist-sat-Vpeak-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",n.transpose([vcirB[:-1], vcirB[1:],  nnV.sum(axis=0)]) )

    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist2d-sat-Mpeak-Accrate-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",dataMC.sum(axis=0))

    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist2d-sat-Vpeak-Accrate-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",dataVC.sum(axis=0))




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
matplotlib.rcParams['font.size']=14
import matplotlib.pyplot as p
import glob
import sys
from scipy.optimize import curve_fit
import cPickle
from os.path import join
from scipy.optimize import minimize

dir = join("D:","data","MultiDark")

def get_cumulative_function(b0, b1, val, volume, minVx = 1e7, maxVx=1e16, Nmin_Occurence=100.):
    """returns the cumulative function n(>X) 
    :param b0: lowerboundary of the bin
    :param b1: higher boundary of the bin
    :param volume:volume of the box in (Mpc/h)^3
    :param minVx: lower value cut. Selectes data in the bins > minVx 
    :param maxVx: higher value cut. Selects data in the bins < maxVx 
    :param Nmin_Occurence: Numberofcounts in a binmustbe larger than Nmin_Occurence tobe considered
    """
    Nc = n.array([ n.sum(val[ii:]) for ii in range(len(val)) ])
    xData_A = (10**b0+10**b1)/2.
    yData_A = Nc/(volume ) 
    yDataErr_A = Nc**(0.5) / volume 
    boundaries = (xData_A > minVx) & (xData_A < maxVx) & (val > Nmin_Occurence )
    sel = (yDataErr_A!=n.inf) & (yData_A>0)  & (boundaries)
    xData = xData_A[sel]
    yData = yData_A[sel]
    yDataErr = yDataErr_A[sel]*yData_A[sel]
    return xData,yData,yDataErr,volume

def get_differential_function(b0, b1, val, volume, minVx = 1e7, maxVx=1e16, Nmin_Occurence=100.):
    """returns the cumulative function n(>X) 
    :param b0: lowerboundary of the bin
    :param b1: higher boundary of the bin
    :param volume:volume of the box in (Mpc/h)^3
    :param minVx: lower value cut. Selectes data in the bins > minVx 
    :param maxVx: higher value cut. Selects data in the bins < maxVx 
    :param Nmin_Occurence: Numberofcounts in a binmustbe larger than Nmin_Occurence tobe considered
    """
    xData_A = (10**b0+10**b1)/2.
    volume_per_bin = (volume * (10**b1-10**b0)) # dV dbin NOT IN LOG
    yData_A = val / volume_per_bin
    yDataErr_A = val**(0.5) / volume_per_bin
    boundaries = (xData_A > minVx) & (xData_A < maxVx) &(val > Nmin_Occurence )
    sel = (yDataErr_A!=n.inf) & (yData_A>0)  & (boundaries)
    xData = xData_A[sel]
    yData = yData_A[sel]
    yDataErr = yDataErr_A[sel]*yData_A[sel]
    return xData,yData,yDataErr,volume

#saundersFct=lambda v, A, logv0,alpha,sig : 10**A * ( v /10**logv0)**(alpha) * n.e**( -n.log10( 1 + v /10**logv0)/(2*sig**2.))
#schechterFct=lambda v, A, logv0,alpha, sig : 10**A * ( v /10**logv0)**(alpha) * n.e**( - v / 10**logv0 /(2*sig**2.) )
#ps * (10**logl/10**logls)**(a+1) * n.e**(-10**logl/10**logls)
#doublePL=lambda v,A,logv0 ,a,b: 10**A * 10**( (1+a)*( v - 10**logv0) + b )


mf = lambda v, A, v0, alpha, beta : 10**A * (v/10**v0)**beta * n.e**(- (v/10**v0)**alpha )

# limits at z0
limits_04 = [1e10, 5e12]
limits_10 = [5e11, 5e13]
limits_25 = [5e12, 5e14]
limits_40 = [1e13, 5e15]
zmin = 0.
zmax = 5

NDecimal = 4

dir_04 = join(dir,"MD_0.4Gpc")
dir_10 = join(dir,"MD_1Gpc")
dir_25 = join(dir,"MD_2.5Gpc")
dir_40 = join(dir,"MD_4Gpc")

dir_boxes =  n.array([dir_04, dir_10, dir_25, dir_40])
zList_files = n.array([ join(dir_box,"redshift-list.txt") for dir_box in dir_boxes])
qty_limits = n.array([limits_04, limits_10, limits_25, limits_40])
volume_boxes =  n.array([400.**3., 1000**3., 2500**3., 4000.**3.])

property_dir = "M200c-mvir"
type = "hist"
cos = "Central" #"Central" # centrak or satellite ?
qty = "M200c"

print "we consider the ",type,"of",qty,"of", cos
print "in the redshift range",zmin,zmax
print "for the boxes",dir_boxes
#print zList_files
print "within the following limits for each box",qty_limits
print "each box has a volume of",volume_boxes, "Mpc3/h3"

fileName = type + "-"+ cos +"-"+ qty +"-*.dat"

############ MD 0.4 Gpc concatenate all functions ##############

fileList = n.array(glob.glob(join(dir_04, property_dir,fileName)))
fileList.sort()

nSN, aSN = n.loadtxt(zList_files[0], unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))

"""
xData_04,yData_04,yDataErr_04,z_04 = [], [], [], []
for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    #print SMDfile
    b0_04, b1_04, val_04 = n.loadtxt(SMDfile,unpack=True)
    xData_04_ii,yData_04_ii,yDataErr_04_ii,volume_04_ii = get_cumulative_function(b0_04, b1_04, val_04,400.**3.,minVx = limits_04[0], maxVx = limits_04[1])
    #print SMDfile.split('-')[-1][:-4]
    z_04_ii = conversion[float(SMDfile.split('-')[-1][:-4])]
    if z_04_ii<zmax and len(xData_04_ii)>0 :
        xData_04.append(xData_04_ii)
        yData_04.append(yData_04_ii)
        yDataErr_04.append(yDataErr_04_ii)
        z_04.append(z_04_ii*n.ones_like(xData_04_ii))

z_04 = n.hstack((z_04))
xData_04 = n.hstack((xData_04))
yData_04 = n.hstack((yData_04))
yDataErr_04 = n.hstack((yDataErr_04))

n.savetxt(join(dir_04, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_cumulative_MD_0.4Gpc.dat"),n.transpose([xData_04,z_04,yData_04,yDataErr_04]), header = qty+" z N Nerr" )

xData_04,yData_04,yDataErr_04,z_04 = [], [], [], []
for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    #print SMDfile
    b0_04, b1_04, val_04 = n.loadtxt(SMDfile,unpack=True)
    xData_04_ii,yData_04_ii,yDataErr_04_ii,volume_04_ii = get_differential_function(b0_04, b1_04, val_04,400.**3.,minVx = limits_04[0], maxVx = limits_04[1])
    #print SMDfile.split('-')[-1][:-4]
    z_04_ii = conversion[float(SMDfile.split('-')[-1][:-4])]
    #print z_04_ii
    if z_04_ii<zmax  and len(xData_04_ii)>0 :
        xData_04.append(xData_04_ii)
        yData_04.append(yData_04_ii)
        yDataErr_04.append(yDataErr_04_ii)
        z_04.append(z_04_ii*n.ones_like(xData_04_ii))

z_04 = n.hstack((z_04))
xData_04 = n.hstack((xData_04))
yData_04 = n.hstack((yData_04))
yDataErr_04 = n.hstack((yDataErr_04))

n.savetxt(join(dir_04, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_differential_MD_0.4Gpc.dat"),n.transpose([xData_04,z_04,yData_04,yDataErr_04]), header = qty+" z N Nerr" )


############ 1 Gpc ##############
fileList = glob.glob(join(dir_10, property_dir,fileName))
nSN, aSN = n.loadtxt(zList_files[1], unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))

xData_10,yData_10,yDataErr_10,z_10 = [], [], [], []
for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    #print SMDfile
    b0_10, b1_10, val_10 = n.loadtxt(SMDfile,unpack=True)
    xData_10_ii,yData_10_ii,yDataErr_10_ii,volume_10_ii = get_cumulative_function(b0_10, b1_10, val_10,1000.**3.,minVx = limits_10[0], maxVx = limits_10[1])
    z_10_ii = conversion[float(SMDfile.split('-')[-1][:-4])]
    if z_10_ii<zmax and len(xData_10_ii)>0 :
        xData_10.append(xData_10_ii)
        yData_10.append(yData_10_ii)
        yDataErr_10.append(yDataErr_10_ii)
        z_10.append(z_10_ii*n.ones_like(xData_10_ii))

z_10 = n.hstack((z_10))
xData_10 = n.hstack((xData_10))
yData_10 = n.hstack((yData_10))
yDataErr_10 = n.hstack((yDataErr_10))

n.savetxt(join(dir_10, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_cumulative_MD_1Gpc"+".dat"),n.transpose([xData_10,z_10,yData_10,yDataErr_10]), header = qty+" z N Nerr")

xData_10,yData_10,yDataErr_10,z_10 = [], [], [], []
for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    #print SMDfile
    b0_10, b1_10, val_10 = n.loadtxt(SMDfile,unpack=True)
    xData_10_ii,yData_10_ii,yDataErr_10_ii,volume_10_ii = get_differential_function(b0_10, b1_10, val_10,1000.**3.,minVx = limits_10[0], maxVx = limits_10[1])
    z_10_ii = conversion[float(SMDfile.split('-')[-1][:-4])]
    if z_10_ii<zmax and len(xData_10_ii)>0 :
        xData_10.append(xData_10_ii)
        yData_10.append(yData_10_ii)
        yDataErr_10.append(yDataErr_10_ii)
        z_10.append(z_10_ii*n.ones_like(xData_10_ii))

z_10 = n.hstack((z_10))
xData_10 = n.hstack((xData_10))
yData_10 = n.hstack((yData_10))
yDataErr_10 = n.hstack((yDataErr_10))

n.savetxt(join(dir_10, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_differential_MD_1Gpc"+".dat"),n.transpose([xData_10,z_10,yData_10,yDataErr_10]), header = qty+" z N Nerr")

############ 2.5 Gpc ##############

fileList = glob.glob(join(dir_25, property_dir,fileName))
nSN, aSN = n.loadtxt(zList_files[2], unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))

xData_25,yData_25,yDataErr_25,z_25 = [], [], [], []
for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    #print SMDfile
    b0_25, b1_25, val_25 = n.loadtxt(SMDfile,unpack=True)
    xData_25_ii,yData_25_ii,yDataErr_25_ii,volume_25_ii = get_cumulative_function(b0_25, b1_25, val_25,2500.**3.,minVx = limits_25[0], maxVx = limits_25[1])
    z_25_ii = conversion[float(SMDfile.split('-')[-1][:-4])]
    if z_25_ii<zmax and len(xData_25_ii):
        xData_25.append(xData_25_ii)
        yData_25.append(yData_25_ii)
        yDataErr_25.append(yDataErr_25_ii)
        z_25.append(z_25_ii*n.ones_like(xData_25_ii))

z_25 = n.hstack((z_25))
xData_25 = n.hstack((xData_25))
yData_25 = n.hstack((yData_25))
yDataErr_25 = n.hstack((yDataErr_25))

n.savetxt(join(dir_25, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_cumulative_MD_2.5Gpc"+".dat"),n.transpose([xData_25,z_25,yData_25,yDataErr_25]), header = qty+" z N Nerr")

xData_25,yData_25,yDataErr_25,z_25 = [], [], [], []
for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    #print SMDfile
    b0_25, b1_25, val_25 = n.loadtxt(SMDfile,unpack=True)
    xData_25_ii,yData_25_ii,yDataErr_25_ii,volume_25_ii = get_differential_function(b0_25, b1_25, val_25,2500.**3.,minVx = limits_25[0], maxVx = limits_25[1])
    z_25_ii = conversion[float(SMDfile.split('-')[-1][:-4])]
    if z_25_ii<zmax and len(xData_25_ii):
        xData_25.append(xData_25_ii)
        yData_25.append(yData_25_ii)
        yDataErr_25.append(yDataErr_25_ii)
        z_25.append(z_25_ii*n.ones_like(xData_25_ii))

z_25 = n.hstack((z_25))
xData_25 = n.hstack((xData_25))
yData_25 = n.hstack((yData_25))
yDataErr_25 = n.hstack((yDataErr_25))

n.savetxt(join(dir_25, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_differential_MD_2.5Gpc"+".dat"),n.transpose([xData_25,z_25,yData_25,yDataErr_25]), header = qty+" z N Nerr")

############ 4 Gpc ##############

fileList = glob.glob(join(dir_40, property_dir,fileName))
nSN, aSN = n.loadtxt(zList_files[3], unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))

xData_40,yData_40,yDataErr_40,z_40 = [], [], [], []
for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    #print SMDfile
    b0_40, b1_40, val_40 = n.loadtxt(SMDfile,unpack=True)
    xData_40_ii,yData_40_ii,yDataErr_40_ii,volume_40_ii = get_cumulative_function(b0_40, b1_40, val_40,4000.**3.,minVx = limits_40[0], maxVx = limits_40[1])
    z_40_ii = conversion[float(SMDfile.split('-')[-1][:-4])]
    if z_40_ii<zmax and len(xData_40_ii):
        xData_40.append(xData_40_ii)
        yData_40.append(yData_40_ii)
        yDataErr_40.append(yDataErr_40_ii)
        z_40.append(z_40_ii*n.ones_like(xData_40_ii))

z_40 = n.hstack((z_40))
xData_40 = n.hstack((xData_40))
yData_40 = n.hstack((yData_40))
yDataErr_40 = n.hstack((yDataErr_40))

n.savetxt(join(dir_40, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_cumulative_MD_4Gpc"+".dat"),n.transpose([xData_40,z_40,yData_40,yDataErr_40]), header = qty+" z N Nerr")

xData_40,yData_40,yDataErr_40,z_40 = [], [], [], []
for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    #print SMDfile
    b0_40, b1_40, val_40 = n.loadtxt(SMDfile,unpack=True)
    xData_40_ii,yData_40_ii,yDataErr_40_ii,volume_40_ii = get_differential_function(b0_40, b1_40, val_40,4000.**3.,minVx = limits_40[0], maxVx = limits_40[1])
    z_40_ii = conversion[float(SMDfile.split('-')[-1][:-4])]
    if z_40_ii<zmax and len(xData_40_ii):
        xData_40.append(xData_40_ii)
        yData_40.append(yData_40_ii)
        yDataErr_40.append(yDataErr_40_ii)
        z_40.append(z_40_ii*n.ones_like(xData_40_ii))

z_40 = n.hstack((z_40))
xData_40 = n.hstack((xData_40))
yData_40 = n.hstack((yData_40))
yDataErr_40 = n.hstack((yDataErr_40))

n.savetxt(join(dir_40, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_differential_MD_4Gpc.dat"),n.transpose([xData_40,z_40,yData_40,yDataErr_40]), header = qty+" z N Nerr")


###########################################################################


property_dir = "M200c-mvir"
type = "hist"
cos = "Satellite" # centrak or satellite ?
qty = "M200c"

print "we consider the ",type,"of",qty,"of", cos
print "in the redshift range",zmin,zmax
print "for the boxes",dir_boxes
#print zList_files
print "within the following limits for each box",qty_limits
print "each box has a volume of",volume_boxes, "Mpc3/h3"

fileName = type + "-"+ cos +"-"+ qty +"-*.dat"

############ MD 0.4 Gpc concatenate all functions ##############

fileList = n.array(glob.glob(join(dir_04, property_dir,fileName)))
fileList.sort()

nSN, aSN = n.loadtxt(zList_files[0], unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))

xData_04,yData_04,yDataErr_04,z_04 = [], [], [], []
for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    #print SMDfile
    b0_04, b1_04, val_04 = n.loadtxt(SMDfile,unpack=True)
    xData_04_ii,yData_04_ii,yDataErr_04_ii,volume_04_ii = get_cumulative_function(b0_04, b1_04, val_04,400.**3.,minVx = limits_04[0], maxVx = limits_04[1])
    #print SMDfile.split('-')[-1][:-4]
    z_04_ii = conversion[float(SMDfile.split('-')[-1][:-4])]
    if z_04_ii<zmax and len(xData_04_ii)>0 :
        xData_04.append(xData_04_ii)
        yData_04.append(yData_04_ii)
        yDataErr_04.append(yDataErr_04_ii)
        z_04.append(z_04_ii*n.ones_like(xData_04_ii))

z_04 = n.hstack((z_04))
xData_04 = n.hstack((xData_04))
yData_04 = n.hstack((yData_04))
yDataErr_04 = n.hstack((yDataErr_04))

n.savetxt(join(dir_04, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_cumulative_MD_0.4Gpc.dat"),n.transpose([xData_04,z_04,yData_04,yDataErr_04]), header = qty+" z N Nerr" )

xData_04,yData_04,yDataErr_04,z_04 = [], [], [], []
for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    #print SMDfile
    b0_04, b1_04, val_04 = n.loadtxt(SMDfile,unpack=True)
    xData_04_ii,yData_04_ii,yDataErr_04_ii,volume_04_ii = get_differential_function(b0_04, b1_04, val_04,400.**3.,minVx = limits_04[0], maxVx = limits_04[1])
    #print SMDfile.split('-')[-1][:-4]
    z_04_ii = conversion[float(SMDfile.split('-')[-1][:-4])]
    #print z_04_ii
    if z_04_ii<zmax  and len(xData_04_ii)>0 :
        xData_04.append(xData_04_ii)
        yData_04.append(yData_04_ii)
        yDataErr_04.append(yDataErr_04_ii)
        z_04.append(z_04_ii*n.ones_like(xData_04_ii))

z_04 = n.hstack((z_04))
xData_04 = n.hstack((xData_04))
yData_04 = n.hstack((yData_04))
yDataErr_04 = n.hstack((yDataErr_04))

n.savetxt(join(dir_04, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_differential_MD_0.4Gpc.dat"),n.transpose([xData_04,z_04,yData_04,yDataErr_04]), header = qty+" z N Nerr" )


############ 1 Gpc ##############
fileList = glob.glob(join(dir_10, property_dir,fileName))
nSN, aSN = n.loadtxt(zList_files[1], unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))

xData_10,yData_10,yDataErr_10,z_10 = [], [], [], []
for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    #print SMDfile
    b0_10, b1_10, val_10 = n.loadtxt(SMDfile,unpack=True)
    xData_10_ii,yData_10_ii,yDataErr_10_ii,volume_10_ii = get_cumulative_function(b0_10, b1_10, val_10,1000.**3.,minVx = limits_10[0], maxVx = limits_10[1])
    z_10_ii = conversion[float(SMDfile.split('-')[-1][:-4])]
    if z_10_ii<zmax and len(xData_10_ii)>0 :
        xData_10.append(xData_10_ii)
        yData_10.append(yData_10_ii)
        yDataErr_10.append(yDataErr_10_ii)
        z_10.append(z_10_ii*n.ones_like(xData_10_ii))

z_10 = n.hstack((z_10))
xData_10 = n.hstack((xData_10))
yData_10 = n.hstack((yData_10))
yDataErr_10 = n.hstack((yDataErr_10))

n.savetxt(join(dir_10, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_cumulative_MD_1Gpc"+".dat"),n.transpose([xData_10,z_10,yData_10,yDataErr_10]), header = qty+" z N Nerr")

xData_10,yData_10,yDataErr_10,z_10 = [], [], [], []
for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    #print SMDfile
    b0_10, b1_10, val_10 = n.loadtxt(SMDfile,unpack=True)
    xData_10_ii,yData_10_ii,yDataErr_10_ii,volume_10_ii = get_differential_function(b0_10, b1_10, val_10,1000.**3.,minVx = limits_10[0], maxVx = limits_10[1])
    z_10_ii = conversion[float(SMDfile.split('-')[-1][:-4])]
    if z_10_ii<zmax and len(xData_10_ii)>0 :
        xData_10.append(xData_10_ii)
        yData_10.append(yData_10_ii)
        yDataErr_10.append(yDataErr_10_ii)
        z_10.append(z_10_ii*n.ones_like(xData_10_ii))

z_10 = n.hstack((z_10))
xData_10 = n.hstack((xData_10))
yData_10 = n.hstack((yData_10))
yDataErr_10 = n.hstack((yDataErr_10))

n.savetxt(join(dir_10, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_differential_MD_1Gpc"+".dat"),n.transpose([xData_10,z_10,yData_10,yDataErr_10]), header = qty+" z N Nerr")

############ 2.5 Gpc ##############

fileList = glob.glob(join(dir_25, property_dir,fileName))
nSN, aSN = n.loadtxt(zList_files[2], unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))

xData_25,yData_25,yDataErr_25,z_25 = [], [], [], []
for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    #print SMDfile
    b0_25, b1_25, val_25 = n.loadtxt(SMDfile,unpack=True)
    xData_25_ii,yData_25_ii,yDataErr_25_ii,volume_25_ii = get_cumulative_function(b0_25, b1_25, val_25,2500.**3.,minVx = limits_25[0], maxVx = limits_25[1])
    z_25_ii = conversion[float(SMDfile.split('-')[-1][:-4])]
    if z_25_ii<zmax and len(xData_25_ii):
        xData_25.append(xData_25_ii)
        yData_25.append(yData_25_ii)
        yDataErr_25.append(yDataErr_25_ii)
        z_25.append(z_25_ii*n.ones_like(xData_25_ii))

z_25 = n.hstack((z_25))
xData_25 = n.hstack((xData_25))
yData_25 = n.hstack((yData_25))
yDataErr_25 = n.hstack((yDataErr_25))

n.savetxt(join(dir_25, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_cumulative_MD_2.5Gpc"+".dat"),n.transpose([xData_25,z_25,yData_25,yDataErr_25]), header = qty+" z N Nerr")

xData_25,yData_25,yDataErr_25,z_25 = [], [], [], []
for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    #print SMDfile
    b0_25, b1_25, val_25 = n.loadtxt(SMDfile,unpack=True)
    xData_25_ii,yData_25_ii,yDataErr_25_ii,volume_25_ii = get_differential_function(b0_25, b1_25, val_25,2500.**3.,minVx = limits_25[0], maxVx = limits_25[1])
    z_25_ii = conversion[float(SMDfile.split('-')[-1][:-4])]
    if z_25_ii<zmax and len(xData_25_ii):
        xData_25.append(xData_25_ii)
        yData_25.append(yData_25_ii)
        yDataErr_25.append(yDataErr_25_ii)
        z_25.append(z_25_ii*n.ones_like(xData_25_ii))

z_25 = n.hstack((z_25))
xData_25 = n.hstack((xData_25))
yData_25 = n.hstack((yData_25))
yDataErr_25 = n.hstack((yDataErr_25))

n.savetxt(join(dir_25, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_differential_MD_2.5Gpc"+".dat"),n.transpose([xData_25,z_25,yData_25,yDataErr_25]), header = qty+" z N Nerr")

############ 4 Gpc ##############

fileList = glob.glob(join(dir_40, property_dir,fileName))
nSN, aSN = n.loadtxt(zList_files[3], unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))

xData_40,yData_40,yDataErr_40,z_40 = [], [], [], []
for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    #print SMDfile
    b0_40, b1_40, val_40 = n.loadtxt(SMDfile,unpack=True)
    xData_40_ii,yData_40_ii,yDataErr_40_ii,volume_40_ii = get_cumulative_function(b0_40, b1_40, val_40,4000.**3.,minVx = limits_40[0], maxVx = limits_40[1])
    z_40_ii = conversion[float(SMDfile.split('-')[-1][:-4])]
    if z_40_ii<zmax and len(xData_40_ii):
        xData_40.append(xData_40_ii)
        yData_40.append(yData_40_ii)
        yDataErr_40.append(yDataErr_40_ii)
        z_40.append(z_40_ii*n.ones_like(xData_40_ii))

z_40 = n.hstack((z_40))
xData_40 = n.hstack((xData_40))
yData_40 = n.hstack((yData_40))
yDataErr_40 = n.hstack((yDataErr_40))

n.savetxt(join(dir_40, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_cumulative_MD_4Gpc"+".dat"),n.transpose([xData_40,z_40,yData_40,yDataErr_40]), header = qty+" z N Nerr")

xData_40,yData_40,yDataErr_40,z_40 = [], [], [], []
for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    #print SMDfile
    b0_40, b1_40, val_40 = n.loadtxt(SMDfile,unpack=True)
    xData_40_ii,yData_40_ii,yDataErr_40_ii,volume_40_ii = get_differential_function(b0_40, b1_40, val_40,4000.**3.,minVx = limits_40[0], maxVx = limits_40[1])
    z_40_ii = conversion[float(SMDfile.split('-')[-1][:-4])]
    if z_40_ii<zmax and len(xData_40_ii):
        xData_40.append(xData_40_ii)
        yData_40.append(yData_40_ii)
        yDataErr_40.append(yDataErr_40_ii)
        z_40.append(z_40_ii*n.ones_like(xData_40_ii))

z_40 = n.hstack((z_40))
xData_40 = n.hstack((xData_40))
yData_40 = n.hstack((yData_40))
yDataErr_40 = n.hstack((yDataErr_40))

n.savetxt(join(dir_40, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_differential_MD_4Gpc.dat"),n.transpose([xData_40,z_40,yData_40,yDataErr_40]), header = qty+" z N Nerr")

"""


xData_04,z_04,yData_04,yDataErr_04 = n.loadtxt(join(dir_04, property_dir,"hist-Central-M200c_ALL_cumulative_MD_0.4Gpc.dat"),unpack=True)
xData_10,z_10,yData_10,yDataErr_10 = n.loadtxt(join(dir_10, property_dir,"hist-Central-M200c_ALL_cumulative_MD_1Gpc.dat"),unpack=True)
xData_25,z_25,yData_25,yDataErr_25 = n.loadtxt(join(dir_25, property_dir,"hist-Central-M200c_ALL_cumulative_MD_2.5Gpc.dat"),unpack=True)
xData_40,z_40,yData_40,yDataErr_40 = n.loadtxt(join(dir_40, property_dir,"hist-Central-M200c_ALL_cumulative_MD_4Gpc.dat"),unpack=True)

s_04 = (z_04 >= zmin) & (z_04 <= zmax)
s_10 = (z_10 >= zmin) & (z_10 <= zmax)
s_25 = (z_25 >= zmin) & (z_25 <= zmax)
s_40 = (z_40 >= zmin) & (z_40 <= zmax)

redshift = n.hstack(( z_04[s_04], z_10[s_10], z_25[s_25], z_40[s_40]))
print "all redshifts available:", set(redshift)
M200c = n.log10(n.hstack(( xData_04[s_04], xData_10[s_10], xData_25[s_25], xData_40[s_40])))
print "min and max masses available:", n.min(M200c), n.max(M200c)
yData = n.log10(n.hstack(( yData_04[s_04], yData_10[s_10], yData_25[s_25], yData_40[s_40])))
print "min and max Y available:", n.min(yData), n.max(yData)
yDataErr = abs(n.hstack(( yDataErr_04[s_04], yDataErr_10[s_10], yDataErr_25[s_25], yDataErr_40[s_40])) / yData)
print "min and max Y error available:", n.min(yDataErr), n.max(yDataErr)


fig = p.figure(1,(5,5))
p.axes([0.19,0.2,0.75,0.75])
sc1=p.scatter(M200c,yData, s=n.ones_like(yData)*3, c=redshift, marker='o',label="MD data", rasterized=True)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'log $M_{200c}$ [h$^{-1}$ M$_{sun}$]')
p.ylim((-9,0))
p.xlim((9,16))
p.ylabel(r'log N$_{cen}$(M>M$_{200c}$) [ h$^3$ Mpc$^{-3}$ ]')
p.grid()
p.show()

xData_04,z_04,yData_04,yDataErr_04 = n.loadtxt(join(dir_04, property_dir,"hist-Central-M200c_ALL_differential_MD_0.4Gpc.dat"),unpack=True)
xData_10,z_10,yData_10,yDataErr_10 = n.loadtxt(join(dir_10, property_dir,"hist-Central-M200c_ALL_differential_MD_1Gpc.dat"),unpack=True)
xData_25,z_25,yData_25,yDataErr_25 = n.loadtxt(join(dir_25, property_dir,"hist-Central-M200c_ALL_differential_MD_2.5Gpc.dat"),unpack=True)
xData_40,z_40,yData_40,yDataErr_40 = n.loadtxt(join(dir_40, property_dir,"hist-Central-M200c_ALL_differential_MD_4Gpc.dat"),unpack=True)

s_04 = (z_04 >= zmin) & (z_04 <= zmax)
s_10 = (z_10 >= zmin) & (z_10 <= zmax)
s_25 = (z_25 >= zmin) & (z_25 <= zmax)
s_40 = (z_40 >= zmin) & (z_40 <= zmax)

redshift = n.hstack(( z_04[s_04], z_10[s_10], z_25[s_25], z_40[s_40]))
print "all redshifts available:", set(redshift)
M200c = n.log10(n.hstack(( xData_04[s_04], xData_10[s_10], xData_25[s_25], xData_40[s_40])))
print "min and max masses available:", n.min(M200c), n.max(M200c)
yData = n.hstack(( yData_04[s_04], yData_10[s_10], yData_25[s_25], yData_40[s_40]))
print "min and max Y available:", n.min(yData), n.max(yData)
yDataErr = n.hstack(( yDataErr_04[s_04], yDataErr_10[s_10], yDataErr_25[s_25], yDataErr_40[s_40])) 
print "min and max Y error available:", n.min(yDataErr), n.max(yDataErr)

hh = aa.H(redshift) /(100 * uu.km / (uu.Mpc* uu.s)) / aa.h
# aa.Om0*
# aa.Om(redshift) *
m2rhom = 10**(2*M200c)/( aa.critical_density(redshift).to(uu.Msun*uu.Mpc**(-3)) * hh**-2)

fig = p.figure(1,(5,5))
p.axes([0.19,0.2,0.75,0.75])
sc1=p.scatter(M200c,m2rhom*yData, s=n.ones_like(yData)*3, c=redshift, marker='o',label="MD data", rasterized=True)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'$\rm \log M_{200c} [h^{-1} M_{sun}$]')
p.ylim((5e-5,0.1))
p.xlim((9,16))
p.yscale('log')
p.ylabel(r'$\rm (M_{200c}^2 /\rho_{cr}(z)) \, dN_{Cen}/dM_{200c}$') #\; [ h^4 Mpc^{-3} M_{sun}^{-1}$]")
p.grid()
p.show()


xData_04,z_04,yData_04,yDataErr_04 = n.loadtxt(join(dir_04, property_dir,"hist-Satellite-M200c_ALL_cumulative_MD_0.4Gpc.dat"),unpack=True)
xData_10,z_10,yData_10,yDataErr_10 = n.loadtxt(join(dir_10, property_dir,"hist-Satellite-M200c_ALL_cumulative_MD_1Gpc.dat"),unpack=True)
#xData_25,z_25,yData_25,yDataErr_25 = n.loadtxt(join(dir_25, property_dir,"hist-Satellite-M200c_ALL_cumulative_MD_2.5Gpc.dat"),unpack=True)
#xData_40,z_40,yData_40,yDataErr_40 = n.loadtxt(join(dir_40, property_dir,"hist-Satellite-M200c_ALL_cumulative_MD_4Gpc.dat"),unpack=True)

s_04 = (z_04 >= zmin) & (z_04 <= zmax)
s_10 = (z_10 >= zmin) & (z_10 <= zmax)
#s_25 = (z_25 >= zmin) & (z_25 <= zmax)
#s_40 = (z_40 >= zmin) & (z_40 <= zmax)

redshift = n.hstack(( z_04[s_04], z_10[s_10])) #, z_25[s_25], z_40[s_40]))
print "all redshifts available:", set(redshift)
M200c = n.log10(n.hstack(( xData_04[s_04], xData_10[s_10] )) )#, xData_25[s_25], xData_40[s_40])))
print "min and max masses available:", n.min(M200c), n.max(M200c)
yData = n.log10(n.hstack(( yData_04[s_04], yData_10[s_10])) )#, yData_25[s_25], yData_40[s_40])))
print "min and max Y available:", n.min(yData), n.max(yData)
yDataErr = abs(n.hstack(( yDataErr_04[s_04], yDataErr_10[s_10])) )#, yDataErr_25[s_25], yDataErr_40[s_40])) / yData)
print "min and max Y error available:", n.min(yDataErr), n.max(yDataErr)

fig = p.figure(1,(5,5))
p.axes([0.19,0.2,0.75,0.75])
sc1=p.scatter(M200c,yData, s=n.ones_like(yData)*3, c=redshift, marker='o',label="MD data", rasterized=True)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'log $M_{200c}$ [h$^{-1}$ M$_{sun}$]')
p.ylim((-9,0))
p.xlim((9,16))
p.ylabel(r'log N$_{sat}$(M>M$_{200c}$) [ h$^3$ Mpc$^{-3}$ ]')
p.grid()
p.show()

xData_04,z_04,yData_04,yDataErr_04 = n.loadtxt(join(dir_04, property_dir,"hist-Satellite-M200c_ALL_differential_MD_0.4Gpc.dat"),unpack=True)
xData_10,z_10,yData_10,yDataErr_10 = n.loadtxt(join(dir_10, property_dir,"hist-Satellite-M200c_ALL_differential_MD_1Gpc.dat"),unpack=True)
xData_25,z_25,yData_25,yDataErr_25 = n.loadtxt(join(dir_25, property_dir,"hist-Satellite-M200c_ALL_differential_MD_2.5Gpc.dat"),unpack=True)
xData_40,z_40,yData_40,yDataErr_40 = n.loadtxt(join(dir_40, property_dir,"hist-Satellite-M200c_ALL_differential_MD_4Gpc.dat"),unpack=True)

s_04 = (z_04 >= zmin) & (z_04 <= zmax)
s_10 = (z_10 >= zmin) & (z_10 <= zmax)
s_25 = (z_25 >= zmin) & (z_25 <= zmax)
s_40 = (z_40 >= zmin) & (z_40 <= zmax)

redshift = n.hstack(( z_04[s_04], z_10[s_10])) #, z_25[s_25], z_40[s_40]))
print "all redshifts available:", set(redshift)
M200c = n.log10(n.hstack(( xData_04[s_04], xData_10[s_10])) )#, xData_25[s_25], xData_40[s_40])))
print "min and max masses available:", n.min(M200c), n.max(M200c)
yData = n.hstack(( yData_04[s_04], yData_10[s_10])) #, yData_25[s_25], yData_40[s_40]))
print "min and max Y available:", n.min(yData), n.max(yData)
yDataErr = n.hstack(( yDataErr_04[s_04], yDataErr_10[s_10])) #, yDataErr_25[s_25], yDataErr_40[s_40])) 
print "min and max Y error available:", n.min(yDataErr), n.max(yDataErr)

hh = aa.H(redshift) /(100 * uu.km / (uu.Mpc* uu.s)) / aa.h
# aa.Om0*
# aa.Om(redshift) *
m2rhom = 10**(2*M200c)/( aa.critical_density(redshift).to(uu.Msun*uu.Mpc**(-3)) * hh**-2)

fig = p.figure(1,(5,5))
p.axes([0.19,0.2,0.75,0.75])
sc1=p.scatter(M200c,m2rhom*yData, s=n.ones_like(yData)*3, c=redshift, marker='o',label="MD data", rasterized=True)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'$\rm \log M_{200c} [h^{-1} M_{sun}$]')
p.ylim((5e-5,0.1))
p.xlim((9,16))
p.yscale('log')
p.ylabel(r'$\rm (M_{200c}^2 /\rho_{cr}(z)) \, dN_{sat}/dM_{200c}$') #\; [ h^4 Mpc^{-3} M_{sun}^{-1}$]")
p.grid()
p.show()
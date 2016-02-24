from os.path import join
import astropy.cosmology as co
#co.Planck13
#co.Planck15
co0307 = co.FlatLambdaCDM(H0=70,Om0=0.3)
#co.WMAP7

import numpy as n

"""
Conversion of bibliography points to Planck 2015 and to the format of data to input the fitting routine.

Needs to match the format:
# z_mean z_min z_max L_mean L_min L_max phi_mean phi_min phi_max  Ngalaxy

No logarithmic values

"""

bibPtDir = join("..", "biblioPoints")
dataDir = join("..", "data")

# 
# load hetdex data. Ciardullo et al. 2013. 0.3 0.7 => Planck 15
names=n.array(["Hetdex-points-z015.dat", "Hetdex-points-z026.dat", "Hetdex-points-z038.dat", "Hetdex-points-z050.dat"])
zmins = n.array([0.0, 0.2, 0.325, 0.45])
zmaxs = n.array([0.2,  0.325, 0.45, 0.56])
zmeans = n.array([0.12,  0.26, 0.38, 0.49])

def convertC13_0(name,zmin,zmax,zmean):
	L,philow,phiup=n.loadtxt(join(bibPtDir,name),unpack=True)
	fV = lambda zMean : co0307.comoving_volume(zMean) / co.Planck15.comoving_volume(zMean)
	fL = lambda zMean : co.Planck15.luminosity_distance(zMean) / co0307.luminosity_distance(zMean)
	z_mean = n.ones_like(philow) * zmean
	z_min = n.ones_like(philow) * zmin
	z_max  = n.ones_like(philow) * zmax
	cL = fL(z_mean)
	L_mean = 10**L * cL
	L_min = 10**(L-0.1) * cL
	L_max = 10**(L+0.1) * cL
	cV = fV(z_mean)
	phi_mean = (phiup*0.001+philow*0.001)/2. * cV
	phi_min = philow*0.001 * cV
	phi_max = phiup*0.001* cV
	Ngalaxy = -1.*n.ones_like(philow) 
	head= "z_mean z_min z_max L_mean L_min L_max phi_mean phi_min phi_max  Ngalaxy"
	f=open(join(dataDir , "O2_3728", "O2_3728-"+name[:-4]+"-Planck15.txt"),'w')
	n.savetxt(f, n.transpose([z_mean, z_min, z_max, L_mean, L_min, L_max, phi_mean, phi_min, phi_max, Ngalaxy]), header= head)
	f.close()

	
def convertC13(name,zmin,zmax,zmean):
	L, phi, philow, phiup=n.loadtxt(join(bibPtDir,name),unpack=True)
	fV = lambda zMean : co0307.comoving_volume(zMean) / co.Planck15.comoving_volume(zMean)
	fL = lambda zMean : co.Planck15.luminosity_distance(zMean) / co0307.luminosity_distance(zMean)
	z_mean = n.ones_like(philow) * zmean
	z_min = n.ones_like(philow) * zmin
	z_max  = n.ones_like(philow) * zmax
	cL = fL(z_mean)
	L_mean = 10**L * cL
	L_min = 10**(L-0.1) * cL
	L_max = 10**(L+0.1) * cL
	cV = fV(z_mean)
	phi_mean = phi*0.001* cV
	phi_min = philow*0.001 * cV
	phi_max = phiup*0.001* cV
	Ngalaxy = -1.*n.ones_like(philow) 
	head= "z_mean z_min z_max L_mean L_min L_max phi_mean phi_min phi_max  Ngalaxy"
	f=open(join(dataDir , "O2_3728", "O2_3728-"+name[:-4]+"-Planck15.txt"),'w')
	n.savetxt(f, n.transpose([z_mean, z_min, z_max, L_mean, L_min, L_max, phi_mean, phi_min, phi_max, Ngalaxy]), header= head)
	f.close()
	
for ii, name in enumerate(names):
	print name, ii
	if ii==0:
		convertC13_0(name,zmins[ii],zmaxs[ii],zmeans[ii])
	if ii>0:
		convertC13(name,zmins[ii],zmaxs[ii],zmeans[ii])


# comparat et al. 2015. Planck 13 => Planck 15
#  Lmin Lmax Lmean phi phiErrTotal jkErrShare poissonErrShare weightErrShare Ngalaxy

def convertC15(name):
	fV = lambda zMean : co.Planck13.comoving_volume(zMean) / co.Planck15.comoving_volume(zMean)
	fL = lambda zMean : co.Planck15.luminosity_distance(zMean) / co.Planck13.luminosity_distance(zMean)
	zmin = float(name[:-4].split('-')[2])
	zmax = float(name[:-4].split('-')[4])
	data=n.loadtxt(join(bibPtDir,name),unpack=True)
	z_mean = n.ones_like(data[0]) * (zmin + zmax)/2.
	z_min = n.ones_like(data[0]) * zmin
	z_max  = n.ones_like(data[0]) * zmax
	cL = fL(z_mean)
	L_mean = 10**data[2] * cL
	L_min = 10**data[0] * cL
	L_max = 10**data[1] * cL
	cV = fV(z_mean)
	phi_mean = data[3] * cV
	phi_min = (data[3] - data[4]) * cV
	phi_max = (data[3] + data[4]) * cV
	Ngalaxy =data[-1]
	head= "z_mean z_min z_max L_mean L_min L_max phi_mean phi_min phi_max  Ngalaxy"
	f=open(join(dataDir , "O2_3728", "O2_3728-"+name[:-4]+"-Planck15.txt"),'w')
	n.savetxt(f, n.transpose([z_mean, z_min, z_max, L_mean, L_min, L_max, phi_mean, phi_min, phi_max, Ngalaxy]), header= head)
	f.close()

names = n.array([ "comparat2014-LF-0.100-z-0.24.dat", "comparat2014-LF-0.240-z-0.4.dat", "comparat2014-LF-0.500-z-0.695.dat", "comparat2014-LF-0.695-z-0.88.dat", "comparat2014-LF-0.880-z-1.09.dat", "comparat2014-LF-1.090-z-1.34.dat", "comparat2014-LF-1.340-z-1.65.dat"])
for name in names:
	convertC15(name)



# load SDSS from Gilbank 2010
sdss01=n.loadtxt("../biblioPoints/Sdss-points-z00.dat",unpack=True)


# the DEEP2 points
dp2=n.loadtxt("../biblioPoints/DEEP2-apj299359t1.txt",unpack=True) # 084 1 12 135
zd2b=n.array([0.84,1.0,1.2,1.35])
jd2=n.array([0,1,2,3])
Lmin=39.5
ymin=5e-5

# Drake te al. 2013
=n.loadtxt(join(bibPtDir,"drake-OII-035-low.csv",unpack=True)
=n.loadtxt(join(bibPtDir,"drake-OII-035-up.csv",unpack=True)
=n.loadtxt(join(bibPtDir,"drake-OII-053-low.csv",unpack=True)
=n.loadtxt(join(bibPtDir,"drake-OII-053-up.csv",unpack=True)
=n.loadtxt(join(bibPtDir,"drake-OII-119-low.csv",unpack=True)
=n.loadtxt(join(bibPtDir,"drake-OII-119-up.csv",unpack=True)
=n.loadtxt(join(bibPtDir,"drake-OII-146-low.csv",unpack=True)
=n.loadtxt(join(bibPtDir,"drake-OII-146-up.csv",unpack=True)
=n.loadtxt(join(bibPtDir,"drake-OII-164-low.csv",unpack=True)
=n.loadtxt(join(bibPtDir,"drake-OII-164-up.csv",unpack=True)
=n.loadtxt(join(bibPtDir,"drake-OIII-063-low.csv",unpack=True)
=n.loadtxt(join(bibPtDir,"drake-OIII-063-up.csv",unpack=True)
=n.loadtxt(join(bibPtDir,"drake-OIII-083-low.csv",unpack=True)
=n.loadtxt(join(bibPtDir,"drake-OIII-083-up.csv",unpack=True)
=n.loadtxt(join(bibPtDir,"drake-OIII-099-low.csv",unpack=True)
=n.loadtxt(join(bibPtDir,"drake-OIII-099-up.csv",unpack=True)


# Zhu et al. 2009
join(bibPtDir,"DEEP2-apj299359t1.txt"

# Ly 2007
join(bibPtDir,"ly-OII-09.dat"
join(bibPtDir,"ly-OIII-04.dat"
join(bibPtDir,"ly-OIII-084-low.csv"
join(bibPtDir,"ly-OIII-084-up.csv"

# Sobral 2012
join(bibPtDir,"sobral-147-fit.dat"
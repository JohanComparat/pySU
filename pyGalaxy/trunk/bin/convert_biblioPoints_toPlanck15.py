from os.path import join
import astropy.cosmology as co
co.Planck13
co.Planck15
co.FlatLambdaCDM(H0=70,Om0=0.3)
co.WMAP7

import numpy as n

"""
Needs to match the format:
# z_mean z_min z_max L_mean L_min L_max phi_mean phi_min phi_max  Ngalaxy
No logarithmic values

"""

bibPtDir = join("..", "biblioPoints")
dataDir = join("..", "data")
fV = lambda zMean : co.Planck13.comoving_volume(zMean) / co.Planck15.comoving_volume(zMean)

fL = lambda zMean : co.Planck15.luminosity_distance(zMean) / co.Planck13.luminosity_distance(zMean)

# comparat et al. 2015. Planck 13 => Planck 15
#  Lmin Lmax Lmean phi phiErrTotal jkErrShare poissonErrShare weightErrShare Ngalaxy
names = n.array([ "comparat2014-LF-0.100-z-0.24.dat", "comparat2014-LF-0.240-z-0.4.dat", "comparat2014-LF-0.500-z-0.695.dat", "comparat2014-LF-0.695-z-0.88.dat", "comparat2014-LF-0.880-z-1.09.dat", "comparat2014-LF-1.090-z-1.34.dat", "comparat2014-LF-1.340-z-1.65.dat"])

def convertC15(name):
	zmin = float(name[:-4].split('-')[2])
	zmax = float(name[:-4].split('-')[4])
	data=n.loadtxt(join(bibPtDir,name),unpack=True)
	z_mean = n.ones_like(data[0]) * (0.1+0.24)/2.
	z_min = n.ones_like(data[0]) * 0.1
	z_max  = n.ones_like(data[0]) * 0.24 
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


for name in names:
	convertC15(name)



# load SDSS
sdss01=n.loadtxt("../biblioPoints/Sdss-points-z00.dat",unpack=True)

# load hetdex data.
hd01=n.loadtxt("../biblioPoints/Hetdex-points-z01.dat",unpack=True)
hd02=n.loadtxt("../biblioPoints/Hetdex-points-z026.dat",unpack=True)
hd03=n.loadtxt("../biblioPoints/Hetdex-points-z038.dat",unpack=True)
hd05=n.loadtxt("../biblioPoints/Hetdex-points-z05.dat",unpack=True)
hdRes=n.array([hd01,hd02,hd03,hd05])
zHd=n.array([0.15,0.26,0.38,0.5])

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
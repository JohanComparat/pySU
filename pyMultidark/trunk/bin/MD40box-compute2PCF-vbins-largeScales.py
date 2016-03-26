from MultiDark import *

box = MultiDarkSimulation(Lbox=4000.0 * uu.Mpc, boxDir = "MD_4Gpc")

all=n.array([ 105, 112, 124, 128, 66, 79, 82, 87, 91, 97 ])
for ii in range(len(all)):
	a= all[ii]
	ll = n.array( glob.glob( "/data2/DATA/eBOSS/Multidark-lightcones/MD_4Gpc/snapshots/hlist_"+str(a)+"_Nb*.fits" ) )
	box.compute2PCF(ll, rmax=200, Nmax=2000000)




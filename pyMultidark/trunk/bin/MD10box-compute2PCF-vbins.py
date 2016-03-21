from MultiDark import *

box = MultiDarkSimulation(Lbox=1000.0 * uu.Mpc, boxDir = "MD_1Gpc_new_rockS")

all=n.array([ 1.00000, 0.50320, 0.60080, 0.74980, 0.89510 ])
for ii in range(len(all)):
	a= all[ii]
	ll = n.array( glob.glob( "/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/*"+str(a)+"*.fits" ) )
	box.compute2PCF(ll)

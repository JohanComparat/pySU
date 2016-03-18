from MultiDark import *

box = MultiDarkSimulation(Lbox=1000.0 * uu.Mpc, boxDir = "MD_1Gpc_new_rockS")

all=n.array([ 0.25320, 0.33030, 0.36090, 0.40320, 0.44060, 0.49220 ])
for ii in range(len(all):
	a= all[ii]
	ll = n.array( glob.glob( "/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/*"+str(a)+"*.fits" ) )
	box.compute2PCF(ll)


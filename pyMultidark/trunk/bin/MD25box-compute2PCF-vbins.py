from MultiDark import *

box = MultiDarkSimulation(Lbox=2500.0 * uu.Mpc, boxDir = "MD_2.5Gpc")

all=n.array([ 46, 74 ])
a= all[0]

ll = n.array( glob.glob( "/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/snapshots/hlist_"+str(a)+"_Nb*.fits" ) )
box.compute2PCF(ll)

a= all[1]

ll = n.array( glob.glob( "/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/snapshots/hlist_"+str(a)+"_Nb*.fits" ) )
box.compute2PCF(ll)



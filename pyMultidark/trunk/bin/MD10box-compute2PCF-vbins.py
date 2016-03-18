from MultiDark import *

box = MultiDarkSimulation(Lbox=1000.0 * uu.Mpc, boxDir = "MD_1Gpc_new_rockS",snl =  n.array(glob.glob(join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_2.5Gpc", "snapshots" , "hlist_*.list"))) ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

all=n.array([ 0.25320, 0.33030, 0.36090 ])
a= all[0]

ll = n.array( glob.glob( "/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/*"+str(a)+"*.fits" ) )
box.compute2PCF(ll)


a= all[1]

ll = n.array( glob.glob( "/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/*"+str(a)+"*.fits" ) )
box.compute2PCF(ll)

a= all[2]

ll = n.array( glob.glob( "/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/*"+str(a)+"*.fits" ) )
box.compute2PCF(ll)

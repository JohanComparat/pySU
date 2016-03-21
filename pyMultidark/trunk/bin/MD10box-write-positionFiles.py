from MultiDark import *

snList= n.array(["/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.50320.list"])#, "/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.74980.list", "/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.89510.list", "/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_1.00000.list" ])

box = MultiDarkSimulation(Lbox=1000.0 * uu.Mpc, boxDir = "MD_1Gpc_new_rockS",snl =snList   ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

for ii in n.arange(len(box.snl)):
	box.writePositionCatalog(ii)


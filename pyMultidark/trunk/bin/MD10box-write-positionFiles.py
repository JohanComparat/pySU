from MultiDark import *

snList= n.array([ "/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.403200.list"])

"""
"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_1.00000.list" , 
"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.89510.list",
"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.74980.list", 
"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.60080.list", 
"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.49220.list", 
"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.44060.list", 
"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.40320.list", 
"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.36090.list", 
"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.33030.list", 
"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.25320.list"
]) 
"""

box = MultiDarkSimulation(Lbox=1000.0 * uu.Mpc, boxDir = "MD_1Gpc_new_rockS",snl =snList   ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

for ii in n.arange(len(box.snl)):
	box.writePositionCatalogPM(ii, logmmin = 11)


from MultiDark import *


snList =  n.array(["/data2/DATA/eBOSS/Multidark-lightcones/MD_4Gpc/snapshots/hlist_66.list", 
"/data2/DATA/eBOSS/Multidark-lightcones/MD_4Gpc/snapshots/hlist_79.list", 
"/data2/DATA/eBOSS/Multidark-lightcones/MD_4Gpc/snapshots/hlist_82.list", 
"/data2/DATA/eBOSS/Multidark-lightcones/MD_4Gpc/snapshots/hlist_87.list", 
"/data2/DATA/eBOSS/Multidark-lightcones/MD_4Gpc/snapshots/hlist_91.list", 
"/data2/DATA/eBOSS/Multidark-lightcones/MD_4Gpc/snapshots/hlist_97.list"])

box = MultiDarkSimulation(Lbox=4000.0 * uu.Mpc, boxDir = "MD_4Gpc",snl = snList ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

for ii in n.arange(len(box.snl)):
	box.writePositionCatalogPM(ii, vmin = 100.)



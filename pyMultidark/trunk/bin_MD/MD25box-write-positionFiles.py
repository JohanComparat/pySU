from MultiDark import *
"""
"/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/snapshots/hlist_10.list",
"/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/snapshots/hlist_74.list",  
"/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/snapshots/hlist_7.list",  
"/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/snapshots/hlist_9.list",
"/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/snapshots/hlist_11.list",
"/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/snapshots/hlist_22.list",
"""
snList =  n.array(["/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/snapshots/hlist_46.list",
"/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/snapshots/hlist_80.list" ])
box = MultiDarkSimulation(Lbox=2500.0 * uu.Mpc, boxDir = "MD_2.5Gpc",snl = snList ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

for ii in n.arange(len(box.snl)):
	box.writePositionCatalogPM(ii, vmin = 100)


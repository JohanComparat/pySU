from MultiDark import *

snList=  n.array([
"/home2/jcomparat/eBOSS-LC/Multidark-lightcones/MD_4GpcNW/snapshots/out_1.list",
"/home2/jcomparat/eBOSS-LC/Multidark-lightcones/MD_4GpcNW/snapshots/out_2.list", 
"/home2/jcomparat/eBOSS-LC/Multidark-lightcones/MD_4GpcNW/snapshots/out_3.list", 
"/home2/jcomparat/eBOSS-LC/Multidark-lightcones/MD_4GpcNW/snapshots/out_4.list", 
"/home2/jcomparat/eBOSS-LC/Multidark-lightcones/MD_4GpcNW/snapshots/out_5.list",
"/home2/jcomparat/eBOSS-LC/Multidark-lightcones/MD_4GpcNW/snapshots/out_6.list",
"/home2/jcomparat/eBOSS-LC/Multidark-lightcones/MD_4GpcNW/snapshots/out_7.list", 
"/home2/jcomparat/eBOSS-LC/Multidark-lightcones/MD_4GpcNW/snapshots/out_8.list", 
"/home2/jcomparat/eBOSS-LC/Multidark-lightcones/MD_4GpcNW/snapshots/out_9.list"])

box = MultiDarkSimulation(Lbox=4000.0 * uu.Mpc, boxDir = "MD_4GpcNW",snl = snList ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

for ii in n.arange(len(box.snl)):
	box.writePositionCatalogPM(ii, vmin = 200., mmin=9*10**11)



from MultiDark import *
snList =  n.array(glob.glob(join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_4Gpc", "snapshots" , "hlist_*.list")))
box = MultiDarkSimulation(Lbox=4000.0 * uu.Mpc, boxDir = "MD_4Gpc",snl = snList ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

for ii in n.arange(len(box.snl)):
	box.writePositionCatalogVmaxM200c(ii, vmin = 200.)



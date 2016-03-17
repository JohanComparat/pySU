from MultiDark import *
box = MultiDarkSimulation(Lbox=2500.0 * uu.Mpc,wdir="/data2/DATA/eBOSS/Multidark-lightcones", boxDir = "MD_2.5Gpc",snl =  n.array(glob.glob(join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_2.5Gpc", "snapshots" , "hlist_*.list"))) ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc),Melement = 23593750000.0)

ll = n.array( glob.glob( "/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/snapshots/*.fits" ) )
vbins = 10**n.arange(n.log10(400), n.log10(3000),0.1)
for ii in n.arange(len(ll)):
	incat = ll[ii]
	for jj in range(len(vbins)-1):
		outfile = incat[:-5] + "_" +str(n.round(vbins[jj],2))+ "_" +str(n.round(vbins[jj+1],2)) + "xiR.pkl"
		box.compute2PCF(incat,outfile,vmin=vbins[jj],vmax=vbins[jj+1],rmin=4,rmax=60)



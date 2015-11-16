from MultiDark import *
box = MultiDarkSimulation(Lbox=2500.0 * uu.Mpc,wdir="/data2/DATA/eBOSS/Multidark-lightcones", boxDir = "MD_1Gpc_new_rockS",snl =  n.array(glob.glob(join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_1Gpc_new_rockS", "snapshots", "hlist_?.?????.list"))) ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc),Melement = 23593750000.0)

for ii in n.arange(len(box.snl))[40::2]:
	box.computeSingleDistributionFunction(ii,'Vpeak', 10**n.arange(0,3.5,0.01))

for ii in n.arange(len(box.snl))[40::2]:
	box.combinesSingleDistributionFunction(ii,'Vpeak', 10**n.arange(0,3.5,0.01),type = "Central")
	box.combinesSingleDistributionFunction(ii,'Vpeak', 10**n.arange(0,3.5,0.01),type = "Satellite")


# parameters 
ii = n.arange(len(box.snl))[40::2][-1]
name='Vpeak' 
bins=10**n.arange(0,3.5,0.01)
type = "Central"

# function
output_dir = join(box.wdir,box.boxDir,"properties",name)
nameSnapshot = box.snl[ii].split('/')[-1][:-5]
pklList = n.array(glob.glob(join(output_dir, nameSnapshot + "_" + name +"_"+type+"_*.pkl")))

nnM = n.empty( [len(pklList),len(bins)-1] ) 
for jj in range(len(pklList)):
    f=open(pklList[jj],'r')
    nnMinter = cPickle.load(f)
    nnM[jj] = nnMinter
    f.close()

n.savetxt(join(output_dir,"hist-"+type+"-"+name+"-"+nameSnapshot[6:]+".dat"),n.transpose([bins[:-1], bins[1:], nnM.sum(axis=0)]))



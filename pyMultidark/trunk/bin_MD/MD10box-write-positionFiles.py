from MultiDark import *

#"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.25320.list",  
#"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.36090.list",  
#"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.44060.list",  
#"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.49220.list",  
#"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.60080.list",  

snList= n.array([ 
"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_1.00000.list", 
"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.89510.list",
"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.33030.list",  
"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.40320.list",  
"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.50320.list",  
"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.74980.list"]) 

box = MultiDarkSimulation(Lbox=1000.0 * uu.Mpc, boxDir = "MD_1Gpc_new_rockS",snl =snList   ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

for ii in n.arange(len(box.snl)):
	box.writePositionCatalogPM(ii,  vmin=50, mmin=2*10**10)

	
snList= n.array(["/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.48140.list"]) 

box = MultiDarkSimulation(Lbox=1000.0 * uu.Mpc, boxDir = "MD_1Gpc_new_rockS",snl =snList   ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))
box.columnDict = {'id': 0, 'desc_id': 1, 'mvir': 2, 'vmax': 3, 'vrms': 4, 'rvir': 5, 'rs': 6, 'Np': 7, 'x': 8, 'y': 9, 'z': 10, 'vx': 11, 'vy': 12, 'vz': 13, 'Jx': 14, 'Jy': 15, 'Jz': 16, 'Spin':17, 'Rs_Klypin': 18, 'Mmvir_all': 19, 'M200b': 20, 'M200c': 21, 'M500c': 22, 'M2500c': 23, 'Xoff': 24, 'Voff': 25, 'Spin_Bullock': 26, 'b_to_a': 27, 'c_to_a': 28, 'Ax': 29, 'Ay': 30, 'Az': 31, 'b_to_a_500c': 32, 'c_to_a_500c': 33, 'Ax_500c': 34, 'Ay_500c': 35, 'Az_500c': 36, 'TU': 37, 'M_pe_Behroozi': 38, 'M_pe_Diemer': 39, 'pid': 40}

for ii in n.arange(len(box.snl)):
	box.writePositionCatalogPM(ii,  vmin=50, mmin=2*10**10)


import os 
import scipy.spatial.ckdtree as t
import time
import numpy as n
from os.path import join
import sys
import glob

Lbox = 1040.
dir = join("PM_Nr1300_L1040.000_g2600_Dgt10.0_ap011_Dlt400._N3em4")
files = glob.glob(join(dir, "GalaxiesZ.????.????.dat"))

readGAL = lambda fileName : n.loadtxt(fileName, unpack=True,skiprows=3)
		
def writeQSOcat(file, dmin=0.3, Lbox = 1040.):
	t0 = time.time()
	print file
	x, y, z, vx, vy, vz, df = readGAL(file)
	catalog = n.transpose([x, y, z, vx, vy, vz*0.92, df])
	
	# identifies pairs
	treeD=t.cKDTree(n.transpose([x, y, z]), 1000.0)
	allpairs = treeD.query_ball_point(treeD.data, dmin)
	# removes auto-pairs
	pairID = []
	for ii, el in enumerate(allpairs):
		if len(el)>1 and ii == el[0]:
			pairID.append(el)

	to_rm = []
	ct = 0
	for id2chooseFrom in pairID:
		n.random.shuffle(id2chooseFrom)
		#keep = id2chooseFrom[0]
		reject = id2chooseFrom[1:]
		to_rm.append(reject)
		ct+=1
	if ct>0:
		to_delete = n.hstack((to_rm)) # list of ids rejected
		print "2del", len(to_delete)
		cat2save = n.delete(catalog, to_delete, axis=0)
		n.savetxt( join(dir, "QSO"+os.path.basename(file)[9:]), cat2save, fmt='%10.5f', header=' x y z vx vy vz df ' )
		print "dt=",time.time()-t0
	else:
		print "no del"
		n.savetxt( join(dir, "QSO"+os.path.basename(file)[9:]), catalog, fmt='%10.5f', header=' x y z vx vy vz df ' )
		print "dt=",time.time()-t0

for file in files:
	if os.path.isfile(join(dir, "QSO"+os.path.basename(file)[9:])):
		pass
	else:
		writeQSOcat(file)


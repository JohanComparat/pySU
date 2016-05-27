#  cd pySU/pyMultidark/trunk/bin/fortranfile-0.2.1/

import numpy as n
import os
from os.path import join
from astropy.io import fits
import time
import fortranfile
DFdir = join("/data2", "users", "gustavo", "BigMD", "1Gpc_3840_Planck1_New", "DENSFIELDS")
mockDir = "/data1/DATA/eBOSS/Multidark-box-mocks/v1.0/parts/"

def writeDFMock(dataCat, DFfile, Lbox = 1000.):
	print dataCat, DFfile
	x, y, z, vx, vy ,vz ,Vmax, Vpeak ,Mvir, parent_id, snap_id, kind, z_dis = n.loadtxt(dataCat ,unpack=True)
	path_to_outputCat = dataCat[:-4] + ".DF.fits.gz"
	# opens the DF file
	f = fortranfile.FortranFile(DFfile)
	gridx, gridy, gridz = f.readInts()
	dx = Lbox/gridx
	# convert QSO positions into indexes
	i = (x/dx).astype(int) 
	j = (y/dx).astype(int) 
	k= (z/dx).astype(int) 
	#init the output array :
	delta = n.empty_like(x)
	delta1 = n.empty_like(x)
	delta2 = n.empty_like(x)
	# now loops over k (every line of the file) and assigns delta values.
	for kk in range(gridx):
		sel = (k==kk)
		N = i[sel] + gridx * j[sel] 
		DF = f.readReals()
		delta[sel] = DF[N]
		# distance 1 mean density field in the plane
		sel1 = (sel)&(i>=1)&(i<gridx-1)&(j>=1)&(j<gridx-1)
		N1 = n.transpose([ (i[sel1]-1) + gridx * (j[sel1] -1), (i[sel1]) + gridx * (j[sel1] -1), (i[sel1]-1) + gridx * (j[sel1]), (i[sel1]+1) + gridx * (j[sel1] +1), (i[sel1]+1) + gridx * (j[sel1] ), (i[sel1]) + gridx * (j[sel1] +1), (i[sel1]+1) + gridx * (j[sel1] -1), (i[sel1]-1) + gridx * (j[sel1] +1) ]) 
		delta1[sel1] = n.array([ n.mean(DF[el]) for el in N1 ]) 
		# assign -1 err value to points on the boundary
		border1 = (sel)&(sel1==False)
		delta1[border1] = n.ones_like(delta1[border1])*-1.
		# distance 2 mean density field in the plane
		sel2 = (sel)&(i>=2)&(i<gridx-2)&(j>=2)&(j<gridx-2)
		N2 = n.transpose([ (i[sel2]-2) + gridx * (j[sel2] -2), (i[sel2]-2) + gridx * (j[sel2] -1), (i[sel2]-2) + gridx * (j[sel2] ), (i[sel2]-2) + gridx * (j[sel2] +1), (i[sel2]-2) + gridx * (j[sel2] +2), (i[sel2]-1) + gridx * (j[sel2] + 2), (i[sel2]) + gridx * (j[sel2] +2), (i[sel2]+11) + gridx * (j[sel2] +2), (i[sel2] + 2) + gridx * (j[sel2] +2), (i[sel2] + 2) + gridx * (j[sel2] +1), (i[sel2] + 2) + gridx * (j[sel2] ), (i[sel2] + 2) + gridx * (j[sel2] -1), (i[sel2] + 2) + gridx * (j[sel2] -2), (i[sel2] + 1) + gridx * (j[sel2] -2),  (i[sel2] ) + gridx * (j[sel2] -2),  (i[sel2] - 1) + gridx * (j[sel2] -2) ]) -1 
		delta2[sel2] = n.array([ n.mean(DF[el]) for el in N2 ])
		# assign -1 err value to points on the boundary
		border2 = (sel)&(sel2==False)
		delta2[border2] = n.ones_like(delta2[border2])*-1.


	f.close()

	c0 = fits.Column(name="DF",format='D', array=delta )
	c01 = fits.Column(name="DF_N1",format='D', array=delta1 )
	c02 = fits.Column(name="DF_N2",format='D', array=delta2 )
	c1 = fits.Column(name="x",format='D', array=x )
	c2 = fits.Column(name="y",format='D', array=y )
	c3 = fits.Column(name="z",format='D', array=z )
	c4 = fits.Column(name="vx",format='D', array=vx )
	c5 = fits.Column(name="vy",format='D', array=vy )
	c6 = fits.Column(name="vz",format='D', array=vz )
	c7 = fits.Column(name="Vmax",format='D', array=Vmax )
	c8 = fits.Column(name="Vpeak",format='D', array=Vpeak )
	c9 = fits.Column(name="Mvir",format='D', array=Mvir )
	c10 = fits.Column(name="parent_id",format='D', array=parent_id )
	c11 = fits.Column(name="snap_id",format='D', array=snap_id )
	c12 = fits.Column(name="kind",format='D', array=kind )
	c13 = fits.Column(name="z_dis",format='D', array=z_dis )
	# now writes the catalog
	cols = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c0, c01, c02 ])
	hdu = fits.BinTableHDU.from_columns(cols)
	os.system("rm -rf "+path_to_outputCat)
	hdu.writeto(path_to_outputCat)


DFfile = join(DFdir,"dmdens_cic_104.dat")
#writeDFMock(join( mockDir,"Box_HAM_z0.701838_nbar1.350000e-05_QSO.dat"), DFfile)
writeDFMock(join( mockDir,"Box_HAM_z0.701838_nbar1.359000e-05_QSO_v1.0.dat"), DFfile)
#writeDFMock(join( mockDir,"Box_HAM_z0.701838_nbar2.400000e-04_ELG.dat"), DFfile)
DFfile = join(DFdir,"dmdens_cic_101.dat")
#writeDFMock(join( mockDir,"Box_HAM_z0.818843_nbar1.000000e-04_LRG.dat"), DFfile)
writeDFMock(join( mockDir,"Box_HAM_z0.818843_nbar1.720000e-05_QSO_v1.0.dat"), DFfile)
#writeDFMock(join( mockDir,"Box_HAM_z0.818843_nbar3.200000e-04_ELG.dat"), DFfile)
DFfile = join(DFdir,"dmdens_cic_097.dat")
writeDFMock(join( mockDir,"Box_HAM_z0.987281_nbar1.933000e-05_QSO_v1.0.dat "), DFfile)
#writeDFMock(join( mockDir,"Box_HAM_z0.987281_nbar2.400000e-04_ELG.dat"), DFfile)
DFfile = join(DFdir,"dmdens_cic_087.dat")
writeDFMock(join( mockDir,"Box_HAM_z1.480160_nbar2.040000e-05_QSO_v1.0.dat"), DFfile)



sys.exit()

N = i[sel] + gridx * j[sel]  # * k[sel]

# distance 1 mean density field in the plane
sel1 = (sel)&(i>=1)&(i<gridx-1)&(j>=1)&(j<gridx-1)

N1 = n.transpose([ (i[sel1]-1) + gridx * (j[sel1] -1), (i[sel1]) + gridx * (j[sel1] -1), (i[sel1]-1) + gridx * (j[sel1]), (i[sel1]+1) + gridx * (j[sel1] +1), (i[sel1]+1) + gridx * (j[sel1] ), (i[sel1]) + gridx * (j[sel1] +1), (i[sel1]+1) + gridx * (j[sel1] -1), (i[sel1]-1) + gridx * (j[sel1] +1) ]) 


delta1[sel1] = n.array([ n.mean(DF[el]) for el in N1 ]) 
# assign -1 err value to points on the boundary
border1 = (sel)&(sel1==False)
delta1[border1] = n.ones_like(delta1[border1])*-1.
"""
# distance 2 mean density field in the plane
sel2 = (sel)&(i>4)&(i<gridx-4)&(j>4)&(j<gridx-4)
N2 = n.transpose([ (i[sel2]-2) * (j[sel2] -2), (i[sel2]-2) * (j[sel2] -1), (i[sel2]-2) * (j[sel2] ), (i[sel2]-2) * (j[sel2] +1), (i[sel2]-2) * (j[sel2] +2), (i[sel2]-1) * (j[sel2] + 2), (i[sel2]) * (j[sel2] +2), (i[sel2]+11) * (j[sel2] +2), (i[sel2] + 2) * (j[sel2] +2), (i[sel2] + 2) * (j[sel2] +1), (i[sel2] + 2) * (j[sel2] ), (i[sel2] + 2) * (j[sel2] -1), (i[sel2] + 2) * (j[sel2] -2), (i[sel2] + 1) * (j[sel2] -2),  (i[sel2] ) * (j[sel2] -2),  (i[sel2] - 1) * (j[sel2] -2) ]) -1 
delta2[sel2] = n.array([ n.mean(DF[el]) for el in N2 ])
# assign -1 err value to points on the boundary
border2 = (sel)&(sel2==False)
delta2[border2] = n.ones_like(delta2[border2])*-1.
"""




ok = (delta > 0)
lgDF = n.log10(delta[ok])

0.9830 has a 0 density field value ?

n.median(lgDF), n.mean(lgDF ), n.min(lgDF), n.max(lgDF)


87 0.403200
97 0.503200
101
104
105 0.600800

0.7 (N=104), 0.85 (N=101), 1 (N=97) and 1.5 (N=87)

 - the first value in the file corresponds to the mean of the density field in the first "pixel" that sits on the point (x, y, z) = (dx/2, dx/2, dx/2) in the box, right ? Where x, y, z are the co moving coordinates given in the snapshot.

 - the second value in the file corresponds to the mean of the density field in the first "pixel" that sits on the point (x, y, z) = (dx/2 + dx, dx/2, dx/2) in the box, right ?

 - the 2048th value is the last one along the x axis, (x, y, z) = (dx/2 + 2047*dx, dx/2, dx/2), then it goes to the second value on the y axis, (x, y, z) = (dx/2 , dx/2 + dx, dx/2), and repeats the operation, right ?

 - in the reading routine, it is not clear if it is first iterated on the x or on the y value. It seems clear that the last value iterated on is the z one : is this right ?


The order of the data in the file goes  as you say:   X,-> Y-> Z. 

In the code, once the  matrix is read into memory, 

RHO(I,J,k) corresponds to  the  point    (i)*dx/2 , (j)*dy/2 ,  (z)*dz/2
where dx = dy=dz = Box/2048



# create fits file

data2 =  np.fromfile(f,dtype="float64",count=134217728)
data3 =  np.fromfile(f,dtype="float64",count=134217728)


data1 =  n.fromfile(f,dtype="float64",count=134217728) # 512 cube
data1 =  n.fromfile(f,dtype="float64",count=134217728) # 512 cube
data1 =  n.fromfile(f,dtype="float64",count=134217728) # 512 cube



x, y, z = 0., 0., 0.
path_to_DF = path_to_file
def get_DF_at_XYZ(x, y, z, path_to_DF, Lbox=1000., gridSize = 2048.):
Lbox=1000.
gridSize = 2048.
dL = Lbox/gridSize

imax = x/dL - 0.5 
imin = x/dL - 0.5 - 1 
jmax = y/dL - 0.5 
jmin = y/dL - 0.5 - 1 
kmax = z/dL - 0.5 
kmin = z/dL - 0.5 - 1 
f=open(path_to_DF,'r')
qty = n.empty( (Nratio,len(bins)-1) ) 
data1 =  n.fromfile(f,dtype="float64",count=NperBatch) # 512 cube


I moved the BOX-mocks to /data1/DATA/eBOSS/Multidark-box-mocks/ In the directory parts/ you can find mocks for each population. In box 104 and 101 we have LRG, ELG and QSO, box 97 has QSO and ELG, and box 87 only has QSO.

The columns are:
#x y z vx vy vz Vmax Vpeak Mvir parent_id snap_id, kind, z_dis

snap_id is the position of the object in the original snapshot, kind = 1 (LRG) 2 (ELG) 3 (QSO), and z_dist, the position in the z direction after applying redshift space distortions.


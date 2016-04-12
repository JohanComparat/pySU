import numpy as n
import os
from os.path import join
from astropy.io import fits
import time

dir = join("/data2", "users", "gustavo", "BigMD", "1Gpc_3840_Planck1_New", "DENSFIELDS")
filename = "dmdens_cic_087.dat"

outDir =  join("/data2", "DATA", "eBOSS", "Multidark-lightcones", "MD_1Gpc_new_rockS", "densityField")

# first rewrite the density field 

path_to_file = join(dir,filename)
f=open(path_to_file,'r')
for i in n.arange(1,65,1):
	print i, time.time()
	path_to_outputCat = join(outDir, "dmdens_cic_087_DF_"+str(i).zfill(2)+".fits.gz")
	data1 =  n.fromfile(f,dtype="float64",count=134217728) # 512 cube
	c0 = fits.Column(name="densityField",format='D', array=data1 )
	cols = fits.ColDefs([c0])
	hdu = fits.BinTableHDU.from_columns(cols)
	os.system("rm -rf "+path_to_outputCat)
	hdu.writeto(path_to_outputCat)


# secondly look for the values of interest

x, y, z, vx, vy ,vz ,Vmax, Vpeak ,Mvir, parent_id, snap_id, kind, z_dis = n.loadtxt("/data1/DATA/eBOSS/Multidark-box-mocks/parts/Box_HAM_z1.480160_nbar1.930000e-05_QSO.dat",unpack=True)
path_to_outputCat = "/data1/DATA/eBOSS/Multidark-box-mocks/parts/Box_HAM_z1.480160_nbar1.930000e-05_QSO.DF.try2.fits.gz"
dx = 1000.0/2048

i = (x/dx).astype(int) + 1
j = (y/dx).astype(int) + 1
k= (z/dx).astype(int) + 1
N = i*j*k
N1 = n.transpose([ 
(i-1)*(j)*(k), 
(i-1)*(j-1)*(k),  (i-1)*(j)*(k-1), 
(i-1)*(j+1)*(k), (i-1)*(j)*(k+1), 
(i-1)*(j-1)*(k-1), (i-1)*(j-1)*(k+1),  (i-1)*(j+1)*(k-1), (i-1)*(j+1)*(k+1), 
(i+1)*(j)*(k), 
(i+1)*(j-1)*(k),  (i+1)*(j)*(k-1), 
(i+1)*(j+1)*(k), (i+1)*(j)*(k+1), 
(i+1)*(j-1)*(k-1), (i+1)*(j-1)*(k+1),  (i+1)*(j+1)*(k-1), (i+1)*(j+1)*(k+1), 
(i)*(j)*(k), 
(i)*(j-1)*(k),  (i)*(j)*(k-1), 
(i)*(j+1)*(k), (i)*(j)*(k+1), 
(i)*(j-1)*(k-1), (i)*(j-1)*(k+1),  (i)*(j+1)*(k-1), (i)*(j+1)*(k+1), 
])

delta = n.empty_like(x)
delta1 = n.empty_like(x)

for i in n.arange(1,65,1):
	print i, time.time()
	path_to_DFCat = join(outDir, "dmdens_cic_087_"+str(i).zfill(2)+".fits.gz")
	hd1 = fits.open( path_to_DFCat )
	sel = (N< i* 512**3)  & (N> (i-1) * 512**3) 
	index = N[ sel ]
	indexD = n.arange(len(x))[sel]
	delta[indexD] = hd1[1].data['densityField'][index - (i-1) * 512**3]
	for ii, el in enumerate(N1):
		sel2 = (el< i* 512**3)  & (el> (i-1) * 512**3) 
		if len(sel2.nonzero()[0])==27:
			delta1[ii] = n.mean(hd1[1].data['densityField'][el - (i-1) * 512**3])

c0 = fits.Column(name="DF",format='D', array=delta )
c01 = fits.Column(name="DF_N1",format='D', array=delta1 )
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

cols = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c0, c01 ])
hdu = fits.BinTableHDU.from_columns(cols)
hdu.writeto(path_to_outputCat)



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


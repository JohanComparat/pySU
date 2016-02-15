import numpy as n
import os
from os.path import join
from astropy.io import fits

dir = join("/data2", "users", "gustavo", "BigMD", "1Gpc_3840_Planck1_New", "DENSFIELDS")
filename = "dmdens_cic_087.dat"

outDir =  join("/data2", "DATA", "eBOSS", "Multidark-lightcones", "MD_1Gpc_new_rockS", "densityField")
path_to_outputCat = join(outDir, "dmdens_cic_087_01.fits")

path_to_file = join(dir,filename)

f=open(path_to_file,'r')
 
 # 64 cubes of 512 on the side ?
data1 =  n.fromfile(f,dtype="float64",count=134217728) # 512 cube

n.histogram(n.log10(data1))

c0 = fits.Column(name="densityField",format='D', array=data1 )
c1 = fits.Column(name="ID",format='I', array=n.arange(len(data1)) )

cols = fits.ColDefs([c0])

hdu = fits.BinTableHDU.from_columns(cols)
hdu.writeto(path_to_outputCat)


# create fits file

data2 =  np.fromfile(f,dtype="float64",count=1000000)
data3 =  np.fromfile(f,dtype="float64",count=1000000)


data1 =  n.fromfile(f,dtype="float64",count=134217728) # 512 cube
data1 =  n.fromfile(f,dtype="float64",count=134217728) # 512 cube
data1 =  n.fromfile(f,dtype="float64",count=134217728) # 512 cube


import numpy as np
import os
from os.path import join

dir = join("/data2", "users", "gustavo", "BigMD", "1Gpc_3840_Planck1_New", "DENSFIELDS")
filename = "dmdens_cic_087.dat"

path_to_file = join(dir,filename)

f=open(path_to_file,'r')
 
data1 =  np.fromfile(f,dtype="float64",count=1 000 000)
data2 =  np.fromfile(f,dtype="float64",count=1 000 000)
data3 =  np.fromfile(f,dtype="float64",count=1 000 000)



from astropy.io import fits
import numpy as n
import glob
import os
import sys

vv = sys.argv[1]

#plates = n.array(os.listdir("/home/comparat/SDSS/26/stellarpop"))
top_dir = "/home/comparat/SDSS/"+vv+"/stellarpop"

#NN = n.array([len(os.listdir('../stellarpop/'+pl)) for pl in plates ])
flyall = n.array(glob.glob(top_dir + "/spFlyPl*.fits"))
flyall.sort()
print(flyall)
print(len(flyall))
#pll = n.array([el.split('/')[-2] for el in flyall])

NN = n.array([len(fits.open(pl)[1].data) for pl in flyall ])
print(NN)
NH = n.array([len(fits.open(pl)[1].data.dtype) for pl in flyall ])
print(set(NH))

n.savetxt(top_dir+'/flyAllList', flyall, fmt='%s')
n.savetxt(top_dir+'/flyAllList_N0', flyall[NN>0], fmt='%s')
n.savetxt(top_dir+'/flyAllList_906', flyall[(NN>0)&(NH==906)], fmt='%s')
n.savetxt(top_dir+'/flyAllList_159', flyall[(NN>0)&(NH==159)], fmt='%s')


import numpy as n
import glob
import os

#plates = n.array(os.listdir("/home/comparat/SDSS/26/stellarpop"))
top_dir = "/home/comparat/SDSS/26/stellarpop"

#NN = n.array([len(os.listdir('../stellarpop/'+pl)) for pl in plates ])
flyall = n.array([glob.glob(top_dir + "/spFlyPl*.fits"))
flyall.sort()

#pll = n.array([el.split('/')[-2] for el in flyall])

NN = n.array([len(fit.open(pl)[1].data) for pl in flyall ])
NH = n.array([len(fit.open(pl)[1].data.dtype) for pl in flyall ])
print(set(NH))

n.savetxt('flyAllList', flyall[NN>0], fmt='%s')
n.savetxt('flyAllList_N0', flyall[NN>0], fmt='%s')
n.savetxt('flyAllList_H0', flyall[NH>0], fmt='%s')
n.savetxt('flyAllList_NH', flyall[(NN>0)&(NH>0)], fmt='%s')


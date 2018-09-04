"""
This script interpolate broad-band photometric filters from different telescopes :
 * DECam : u, g, r, i, z
 * CFHT / MEGACAM : u, g, r, i, z
 * SDSS : u, g, r, i, z

It also provides the mean wavelength of each filter.

"""
from os.path import join
import os
import numpy as n
from scipy.interpolate import interp1d
from scipy.integrate import quad

filterDir = join(os.environ['GIT_PYSU'],"galaxy/data","photometricFilterDir")

dt=n.loadtxt(join(filterDir, "decamFilter","decam_u.par"),unpack=True)
filterUdecam=interp1d(n.hstack((2000.,dt[0],12000.)),n.hstack((0.,dt[1],0.)))
toInt=interp1d(dt[0],dt[0]*dt[1])
num=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
toInt=interp1d(dt[0],dt[1])
den=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
lambUdecam=num[0]/den[0]  # 3816.9196881965267

dt=n.loadtxt(join(filterDir, "decamFilter","decam_g.par"),unpack=True)
filterGdecam=interp1d(n.hstack((2000.,dt[0],12000.)),n.hstack((0.,dt[1],0.)))
toInt=interp1d(dt[0],dt[0]*dt[1])
num=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
toInt=interp1d(dt[0],dt[1])
den=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
lambGdecam=num[0]/den[0]  # 4826.595621922894

dt=n.loadtxt(join(filterDir, "decamFilter","decam_r.par"),unpack=True)
filterRdecam=interp1d(n.hstack((2000.,dt[0],12000.)),n.hstack((0.,dt[1],0.)))
toInt=interp1d(dt[0],dt[0]*dt[1])
num=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
toInt=interp1d(dt[0],dt[1])
den=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
lambRdecam=num[0]/den[0]  # 6435.200607130994

dt=n.loadtxt(join(filterDir, "decamFilter","decam_i.par"),unpack=True)
filterIdecam=interp1d(n.hstack((2000.,dt[0],12000.)),n.hstack((0.,dt[1],0.)))
toInt=interp1d(dt[0],dt[0]*dt[1])
num=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
toInt=interp1d(dt[0],dt[1])
den=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
lambIdecam=num[0]/den[0]  # 7825.443132094219

dt=n.loadtxt(join(filterDir, "decamFilter","decam_z.par"),unpack=True)
filterZdecam=interp1d(n.hstack((2000.,dt[0],12000.)),n.hstack((0.,dt[1],0.)))
toInt=interp1d(dt[0],dt[0]*dt[1])
num=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
toInt=interp1d(dt[0],dt[1])
den=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
lambZdecam=num[0]/den[0]  # 9179.697101576936


dt=n.loadtxt(join(filterDir, "cfhtFilter","megacamQE.dat"),unpack=True,usecols=(0,1))
qe=interp1d(n.hstack((2000.,dt[0],12000.)),n.hstack((0.,dt[1]/100.,0.)))

wl,u,g,i,r,z=n.loadtxt(join(filterDir, "cfhtFilter","megacamFilters.dat"),unpack=True)
xx=n.hstack((2000.,wl[n.argsort(wl)]*10.,12000.))
filterUcfht=interp1d(xx,n.hstack((0.,u[n.argsort(wl)]/100.,0.))*qe(xx))
filterGcfht=interp1d(xx,n.hstack((0.,g[n.argsort(wl)]/100.,0.))*qe(xx))
filterRcfht=interp1d(xx,n.hstack((0.,r[n.argsort(wl)]/100.,0.))*qe(xx))
filterIcfht=interp1d(xx,n.hstack((0.,i[n.argsort(wl)]/100.,0.))*qe(xx))
filterZcfht=interp1d(xx,n.hstack((0.,z[n.argsort(wl)]/100.,0.))*qe(xx))

toInt=interp1d(filterUcfht.x,filterUcfht.x*filterUcfht.y)
num=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
toInt=interp1d(filterUcfht.x,filterUcfht.y)
den=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
lambUcfht=num[0]/den[0]  # 3798.5222429284845

toInt=interp1d(filterGcfht.x,filterGcfht.x*filterGcfht.y)
num=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
toInt=interp1d(filterGcfht.x,filterGcfht.y)
den=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
lambGcfht=num[0]/den[0]  # 4863.631644522178

toInt=interp1d(filterRcfht.x,filterRcfht.x*filterRcfht.y)
num=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
toInt=interp1d(filterRcfht.x,filterRcfht.y)
den=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
lambRcfht=num[0]/den[0]  # 6262.08456842134

toInt=interp1d(filterIcfht.x,filterIcfht.x*filterIcfht.y)
num=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
toInt=interp1d(filterIcfht.x,filterIcfht.y)
den=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
lambIcfht=num[0]/den[0]  # 7683.83577564458

toInt=interp1d(filterZcfht.x,filterZcfht.x*filterZcfht.y)
num=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
toInt=interp1d(filterZcfht.x,filterZcfht.y)
den=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
lambZcfht=num[0]/den[0]  # 9030.5554905110


dt=n.loadtxt(join(filterDir, "sdssFilter","up.pb"),unpack=True)
filterUsdss=interp1d(n.hstack((2000.,dt[0],12000.)),n.hstack((0.,dt[1],0.)))
toInt=interp1d(filterUsdss.x,filterUsdss.x*filterUsdss.y)
num=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
toInt=interp1d(filterUsdss.x,filterUsdss.y)
den=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
lambUsdss=num[0]/den[0]  # 6185.194441906053

dt=n.loadtxt(join(filterDir, "sdssFilter","gp.pb"),unpack=True)
filterGsdss=interp1d(n.hstack((2000.,dt[0],12000.)),n.hstack((0.,dt[1],0.)))
toInt=interp1d(filterGsdss.x,filterGsdss.x*filterGsdss.y)
num=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
toInt=interp1d(filterGsdss.x,filterGsdss.y)
den=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
lambGsdss=num[0]/den[0]  # 6185.194441906053

dt=n.loadtxt(join(filterDir, "sdssFilter","rp.pb"),unpack=True)
filterRsdss=interp1d(n.hstack((2000.,dt[0],12000.)),n.hstack((0.,dt[1],0.)))
toInt=interp1d(filterRsdss.x,filterRsdss.x*filterRsdss.y)
num=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
toInt=interp1d(filterRsdss.x,filterRsdss.y)
den=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
lambRsdss=num[0]/den[0]  # 6185.194441906053

dt=n.loadtxt(join(filterDir, "sdssFilter","ip.pb"),unpack=True)
filterIsdss=interp1d(n.hstack((2000.,dt[0],12000.)),n.hstack((0.,dt[1],0.)))
toInt=interp1d(filterIsdss.x,filterIsdss.x*filterIsdss.y)
num=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
toInt=interp1d(filterIsdss.x,filterIsdss.y)
den=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
lambIsdss=num[0]/den[0]  # 6185.194441906053

dt=n.loadtxt(join(filterDir, "sdssFilter","zp.pb"),unpack=True)
filterZsdss=interp1d(n.hstack((2000.,dt[0],12000.)),n.hstack((0.,dt[1],0.)))
toInt=interp1d(filterZsdss.x,filterZsdss.x*filterZsdss.y)
num=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
toInt=interp1d(filterZsdss.x,filterZsdss.y)
den=quad(toInt,toInt.x.min()+10,toInt.x.max()-10,limit=500000)
lambZsdss=num[0]/den[0]  # 6185.194441906053





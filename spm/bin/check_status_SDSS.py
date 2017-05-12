
from os.path import join
import os
import numpy as n
import glob 
import sys 
import time
import astropy.io.fits as fits

env = os.environ['SDSSDR12_DIR']
path_2_file = join(env, "catalogs", "specObj-dr12.fits")
data = fits.open(path_2_file)[1].data

selection = (data['ZWARNING']==0) & (data['CLASS']=="GALAXY") & (data['Z'] > data['Z_ERR']) & (data['Z_ERR']>0) 

dirs = n.array([
"", 
"stellarpop-m11-chabrier-miles", "stellarpop-m11-chabrier-stelib", "stellarpop-m11-chabrier-elodie", 
"stellarpop-m11-kroupa-miles",   "stellarpop-m11-kroupa-stelib",  "stellarpop-m11-kroupa-elodie", 
"stellarpop-m11-salpeter-miles", "stellarpop-m11-salpeter-stelib",   "stellarpop-m11-salpeter-elodie"])  
suffixes = n.array([
".fits", 
"-cha.fits", "-cha.fits", "-cha.fits", 
"-kr.fits", "-kr.fits", "-kr.fits", 
"-ss.fits", "-ss.fits", "-ss.fits"])

done=n.zeros((len(data), len(dirs)))

bds = n.arange(0,len(data),100000)
for jj in n.arange(0,len(data)+100000,100000)[:-1]:
	for ii, el in enumerate(data[bds[jj]:bds[jj+1]]) :
		plate = str(int(el['PLATE'])).zfill(4)
		mjd = str(int(el['MJD']))
		fiber = str(int(el['FIBERID'])).zfill(4)
		flyName = "spFly-"+plate+"-"+mjd+"-"+fiber
		done[ii] = 1.*n.array([ os.path.isfile(join(env, dir, 'stellarpop', plate, flyName + suff)) for dir, suff in zip(dirs, suffixes) ])
		print done[ii]

	outData = n.vstack((data['PLATE'], data['MJD'], data['FIBERID'], done.T)).astype('int')

	n.savetxt(os.environ['DATA_DIR'], 'status', 'status-SDSS-'+str(bds[jj])+'.txt.gz', outData[selection], fmt='%i', header=' plate mjd fiberid summary cm cs ce km ks ke sm ss se ' )


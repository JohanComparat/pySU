#!/bin/bash
#PBS -l walltime=4:00:00
#PBS -o ebossdr14_kroupa.o.$PBS_JOBID
#PBS -e ebossdr14_kroupa.e$PBS_JOBID
#PBS -M comparat@mpe.mpg.de
module load apps/anaconda/2.4.1
module load apps/python/2.7.8/gcc-4.4.7

cd /users/comparat/batch_spFly_perPlate
python createPlSpAll.py

from os.path import join
import os
import numpy as n
import glob 
import sys 
import time
import astropy.io.fits as fits

init_cat = join( os.environ['EBOSSDR14_DIR'], "catalogs", "spAll-v5_10_0.fits")
hdus = fits.open(init_cat)
tbdata = hdus[1].data
plates = n.array(list(set(tbdata['PLATE'] )))
plates.sort()

for plate in plates:
	print plate, time.time()
	mask = tbdata['PLATE'] ==int(plate)
	if len(mask.nonzero()[0])>=1 :
		newtbdata = tbdata[mask]
		hdu = fits.BinTableHDU(data=newtbdata)
		newTab=join( os.environ['EBOSSDR14_DIR'], "catalogs", "perPlate", "sp-"+str(plate).zfill(4)+".fits")
		hdu.writeto(newTab)

		
init_cat = join( os.environ['SDSSDR12_DIR'], "catalogs", "specObj-SDSS-dr12.fits")
hdus = fits.open(init_cat)
tbdata = hdus[1].data
plates = n.array(list(set(tbdata['PLATE'] )))
plates.sort()

for plate in plates:
	print plate, time.time()
	mask = tbdata['PLATE'] ==int(plate)
	if len(mask.nonzero()[0])>=1 :
		newtbdata = tbdata[mask]
		hdu = fits.BinTableHDU(data=newtbdata)
		newTab=join( os.environ['SDSSDR12_DIR'], "catalogs", "perPlate", "sp-"+str(plate).zfill(4)+".fits")
		hdu.writeto(newTab)

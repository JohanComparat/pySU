import glob
import os
from os.path import join
import numpy as n

def writeScript(rootName, plate, env):
	f=open(rootName+".sh",'w')
	f.write("#!/bin/bash \n")
	f.write("#PBS -l walltime=4:00:00 \n")
	f.write("#PBS -o "+plate+".o.$PBS_JOBID \n")
	f.write("#PBS -e "+plate+".e$PBS_JOBID \n")
	f.write("#PBS -M comparat@mpe.mpg.de \n")
	f.write("module load apps/anaconda/2.4.1 \n")
	f.write("module load apps/python/2.7.8/gcc-4.4.7 \n")
	f.write("export PYTHONPATH=$PYTHONPATH:/users/comparat/pySU/galaxy/python/ \n")
	f.write("export PYTHONPATH=$PYTHONPATH:/users/comparat/pySU/spm/python/ \n")
	f.write(" \n")
	f.write("cd /users/comparat/pySU/spm/bin \n")

	specList = n.array(glob.glob(os.path.join(os.environ[env], 'stellarpop-m11-chabrier', 'stellarpop', plate, 'spFly*.fits')))
	data = n.array([os.path.basename(specName).split('-') for specName in specList]) 
	
	for el in data :
		f.write("python combine_model_spectra.py "+el[1]+" "+el[2]+" "+el[3]+" "+env+" \n")
	
	f.write(" \n")
	f.close()


env="EBOSSDR14_DIR"
plates = n.loadtxt( join(os.environ[env], "catalogs", "plateNumberList"), unpack=True, dtype='str')
for plate in plates:
	rootName = join(os.environ['HOME'], "batch_combine", plate)
	writeScript(rootName, plate, env)

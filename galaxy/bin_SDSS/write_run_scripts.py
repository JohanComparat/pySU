import os
from os.path import join
import numpy as n
def writeScript(rootName, plate):
	f=open(rootName+".sh",'w')
	f.write("#!/bin/bash \n")
	f.write("#PBS -l walltime=48:00:00 \n")
	f.write("#PBS -o "+plate+".o.$PBS_JOBID \n")
	f.write("#PBS -e "+plate+".e$PBS_JOBID \n")
	f.write("#PBS -M johan.comparat@gmail.com \n")
	f.write("module load apps/anaconda/2.4.1 \n")
	f.write("module load apps/python/2.7.8/gcc-4.4.7 \n")
	f.write("export PYTHONPATH=$PYTHONPATH:/users/comparat/pySU/galaxy/python/ \\n")
	f.write("export PYTHONPATH=$PYTHONPATH:/users/comparat/pySU/simulations/python/ \n")
	f.write("export PYTHONPATH=$PYTHONPATH:/users/comparat/pySU/multidark/python/ \n")
	f.write("export PYTHONPATH=$PYTHONPATH:/users/comparat/pySU/spm/python/ \n")
	f.write("export PYTHONPATH=$PYTHONPATH:/users/comparat/pySU/targetselection/python/ \n")
	f.write(" \n")
	f.write("cd /users/comparat/pySU/galaxy/bin_SDSS \n")
	f.write("python model_stacked_spectra_SDSS "+plate+" \n")
	f.write(" \n")
	f.close()

plateList = n.loadtxt("plateList", unpack = True)
for plate in plateList:
	rootName = join(os.environ['HOME'], "batchscripts_firefly_model", plate)
	writeScript(rootName, plate)
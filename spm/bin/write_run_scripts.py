import os
plate = "7645"
rootName = join(os.environ['HOME'], "batchscripts_firefly", plate)

f=open(rootName+".sh",'w')

str2write = """
#!/bin/bash 
#PBS -l walltime=48:00:00 
#PBS -o """+plate+""".o.$PBS_JOBID 
#PBS -e """+plate+""".e$PBS_JOBID 
#PBS -M johan.comparat@gmail.com 

# Configure modules and load modules 
module load apps/anaconda/2.4.1
module load apps/python/2.7.8/gcc-4.4.7

export PYTHONPATH=$PYTHONPATH:/users/comparat/pySU/galaxy/python/
export PYTHONPATH=$PYTHONPATH:/users/comparat/pySU/simulations/python/
export PYTHONPATH=$PYTHONPATH:/users/comparat/pySU/multidark/python/
export PYTHONPATH=$PYTHONPATH:/users/comparat/pySU/spm/python/
export PYTHONPATH=$PYTHONPATH:/users/comparat/pySU/targetselection/python/

cd /users/comparat/pySU/spm/bin
python stellarpop_sdss_single """+plate

f.write(str2write)
f.close()

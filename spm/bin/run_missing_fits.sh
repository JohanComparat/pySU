#!/bin/bash 
#PBS -l walltime=90:00:00 
#PBS -o missing.o.$PBS_JOBID 
#PBS -e missing.e.$PBS_JOBID 
#PBS -M johan.comparat@gmail.com 
module load apps/anaconda/2.4.1 
module load apps/python/2.7.8/gcc-4.4.7 
export PYTHONPATH=$PYTHONPATH:/users/comparat/pySU/galaxy/python/
export PYTHONPATH=$PYTHONPATH:/users/comparat/pySU/spm/python/ 
 
cd /users/comparat/pySU/spm/bin 

python stellarpop_sdss_singleSpec_kroupa 274 51913 66

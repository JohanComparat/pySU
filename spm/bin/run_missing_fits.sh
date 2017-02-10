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

python stellarpop_sdss_singleSpec_kroupa 267 51608 50
python stellarpop_sdss_singleSpec_kroupa 267 51608 52
python stellarpop_sdss_singleSpec_kroupa 267 51608 54
python stellarpop_sdss_singleSpec_kroupa 267 51608 55
python stellarpop_sdss_singleSpec_kroupa 267 51608 56
python stellarpop_sdss_singleSpec_kroupa 267 51608 57
python stellarpop_sdss_singleSpec_kroupa 267 51608 58
python stellarpop_sdss_singleSpec_kroupa 267 51608 59
python stellarpop_sdss_singleSpec_kroupa 267 51608 60
python stellarpop_sdss_singleSpec_kroupa 267 51608 61
python stellarpop_sdss_singleSpec_kroupa 267 51608 62
python stellarpop_sdss_singleSpec_kroupa 267 51608 63
python stellarpop_sdss_singleSpec_kroupa 267 51608 64
python stellarpop_sdss_singleSpec_kroupa 267 51608 65
python stellarpop_sdss_singleSpec_kroupa 267 51608 66
python stellarpop_sdss_singleSpec_kroupa 274 51913 49
python stellarpop_sdss_singleSpec_kroupa 274 51913 51
python stellarpop_sdss_singleSpec_kroupa 274 51913 53
python stellarpop_sdss_singleSpec_kroupa 274 51913 54
python stellarpop_sdss_singleSpec_kroupa 274 51913 56
python stellarpop_sdss_singleSpec_kroupa 274 51913 58
python stellarpop_sdss_singleSpec_kroupa 274 51913 59
python stellarpop_sdss_singleSpec_kroupa 274 51913 62
python stellarpop_sdss_singleSpec_kroupa 274 51913 65
python stellarpop_sdss_singleSpec_kroupa 274 51913 66

#!/bin/bash
#SBATCH -o /home2/jcomparat/eBOSS-LC/PM_runs/PM_C/PM_L650_g1612_Dgt120_am100_Dlt240/test_job.%j.%N.out
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=all
#SBATCH --export=NONE
#SBATCH --time=167:00:00
#SBATCH -A 32cores 
#SBATCH --mem=700000
#SBATCH --job-name=PM04

export OMP_NUM_THREADS=32

echo 1 | ../PMP2start.exe 
echo 1000000 | ../PMP2main.exe


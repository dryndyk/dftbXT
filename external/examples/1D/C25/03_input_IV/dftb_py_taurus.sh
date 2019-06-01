#!/bin/bash
#SBATCH -J C25_03
#SBATCH --time=24:00:00
#SBATCH -o out
#SBATCH --ntasks=100
##SBATCH --nodes=13
##SBATCH --ntasks-per-node=8
#SBATCH	--mail-type=all
#SBATCH -A p_nanocarbon
#SBATCH --mem-per-cpu=7500
##SBATCH --mem=60000
#SBATCH -p haswell

echo Starting Program

module load modenv/classic
module load gcc/6.3.0 openmpi/2.1.0-gnu6.3
source /home/dryndyk/nanocarbon/codes/cp2kXT/cp2k/tools/toolchain/install/setup
export OMP_NUM_THREADS=1
export PATH=$PATH:~/bin

python IV.py

echo Finished Program 



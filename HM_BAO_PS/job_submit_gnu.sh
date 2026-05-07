#!/bin/bash
#
#SBATCH -J HM_PS_Calc

#SBATCH -N 4
#SBATCH --ntasks-per-node=16

#SBATCH -p long
##SBATCH -p normal
#SBATCH --mem=100G
#SBATCh --exclusive
#SBATCH --no-requeue 

cd $SLURM_SUBMIT_DIR
module purge

module load compiler/devtoolset/7.3.1
module load compiler/rocm/dtk/22.10.1
module load mpi/hpcx/2.11.0/gcc-7.3.1
module load mathlib/hdf5/1.8.20/gcc


srun --mpi=pmix_v3 ./Power.exe ./HM4T_PS_Halo.tex 3

date


################################
#$ -S /bin/bash
#$ -j y
#$ -cwd
#$ -m n
### must be a multiple of 16
#$ -pe impi_hydra 64
#$ -l h_rt=24:00:00

module load impi

#time mpiexec -perhost 4 -n 32 ./L-VelHalo-MXXL mxxl_dm.tex 54

time mpiexec -perhost 8 -n 32 ./L-VelHalo-MXXL mxxl_dm_field.tex 54
#time mpiexec -perhost 8 -n 32 ./L-VelHalo-MXXL mxxl_halo_field.tex 54
#time mpiexec -perhost 8 -n 32 ./L-VelHalo-MXXL mxxl_subhalo_field.tex 54
date


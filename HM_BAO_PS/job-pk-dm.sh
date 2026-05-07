#$ -S /bin/bash
#$ -j y
#$ -cwd
#$ -m n
### must be a multiple of 12
#$ -pe impi4 132 
### up to 12 hours
#$ -l h_rt=24:00:00

module load impi/4.0.3

mpiexec -np 128 ./Power.exe mxxl_dm.tex 0


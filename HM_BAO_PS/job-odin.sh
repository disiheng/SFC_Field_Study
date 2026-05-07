#$ -S /bin/bash
#$ -j y
#$ -cwd
#$ -m n
### must be a multiple of 12
#$ -pe impi_hydra 128
### up to 12 hours
#$ -l h_rt=23:00:00

module load impi

mpiexec -perhost 16 -np 128 ./Power.exe mill_om015.tex 7 
#mpiexec -perhost 16 -np 256 ./Power.exe mill_om025.tex 7 
#mpiexec -perhost 16 -np 256 ./Power.exe mill_om029.tex 7 
#mpiexec -perhost 16 -np 256 ./Power.exe mill_om05.tex 7 
#mpiexec -perhost 16 -np 256 ./Power.exe mill_om08.tex 7 

#mpiexec -np 128 ./Power.exe mxxl_gals.tex 1 ${SGE_TASK_ID}
#mpiexec -np 128 ./Power.exe mxxl_gals.tex 2 ${SGE_TASK_ID}
#mpiexec -np 128 ./Power.exe mxxl_gals.tex 1 
#mpiexec -np 128 ./Power.exe mill_gals.tex 1 
#piexec -np 128 ./Power.exe mxxl_halos.tex 3
#mpiexec -np 128 ./Power.exe mxxl_halos.tex 4
#mpiexec -np 128 ./Power.exe mxxl_halos.tex 5

#mpiexec -np 32 ./Power.exe mill_om015.tex 0 ${SGE_TASK_ID} 
#mpiexec -np 32 ./Power.exe mill_om025.tex 0 ${SGE_TASK_ID} 
#mpiexec -np 32 ./Power.exe mill_om05.tex 0 ${SGE_TASK_ID} 
#mpiexec -np 32 ./Power.exe mill_om08.tex 0 ${SGE_TASK_ID} 
#mpiexec -perhost 8 -np 32 ./Power.exe mII_dm.tex 0 ${SGE_TASK_ID}
#mpiexec -perhost 8 -np 32 ./Power.exe mill_large_om05.tex 0 ${SGE_TASK_ID} 
#mpiexec -np 32 ./Power.exe mill_om08.tex 0 ${SGE_TASK_ID} 
#mpiexec -np 64 ./Power.exe lbassic_dm.tex 0 


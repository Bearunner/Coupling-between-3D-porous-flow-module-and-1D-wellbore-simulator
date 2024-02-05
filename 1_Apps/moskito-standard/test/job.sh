#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=28
#SBATCH --time=00:05:00
#SBATCH --export=ALL,PATH_TO=/projects/moskito/,EXECUTABLE=moskito-opt,INPUT=2c.i

#module load ${MPI_MODULE}
#module load numlib/petsc/3.8.3-openmpi-3.1-gnu-8.2 # Not needed, if custom PETSc installation used 

cd ${SLURM_SUBMIT_DIR}
#cd ~/projects/moskito
startexe="mpirun --bind-to core --map-by core ${HOME}${PATH_TO}${EXECUTABLE} -i ${SLURM_SUBMIT_DIR}/${INPUT}"

echo $startexe
exec $startexe

exit
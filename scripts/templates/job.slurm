#!/bin/bash -x
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=48
#SBATCH --output=mpi-out.%j
#SBATCH --error=mpi-err.%j
#SBATCH --time=02:00:00


# http://www-01.ibm.com/support/docview.wss?uid=isg3T1024294
# Best practive for LARGE parallel jobs with intelMPI 
# DOES NOT WORK! 
#export I_MPI_HYDRA_BOOTSTRAP=lsf # Set to lsf means to use LSF blaunch
#export I_MPI_HYDRA_BRANCH_COUNT=-1 # The default value is -1 if less than 128 nodes are used. This value also means that there is no hierarchical structure
#export I_MPI_LSF_USE_COLLECTIVE_LAUNCH=1 # Set to 1 to let Intel MPI use blaunch ?z to launch tasks.

# buffer sizes for infiniband; 
# keep small to reduce memory overhead for large parallel jobs
export I_MPI_DAPL_BUFFER_NUM=4
export I_MPI_DAPL_BUFFER_SIZE=128

export OMP_NUM_THREADS=2
export OMP_MAX_THREADS=2 
export KMP_AFFINITY=compact

#IntelMPI process pinning 
export I_MPI_PIN=yes       # force process pinning 
export I_MPI_PIN_MODE=pm   # use pinning by launcher; should collocate CPU and memory; not guaranteed for option 'lib' 
export I_MPI_PIN_DOMAIN=omp 


# ./dns_serial.x
srun ./dns_fbc.x 
#srun ./dns_omp.x

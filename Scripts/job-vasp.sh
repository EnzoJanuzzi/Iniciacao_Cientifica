#!/bin/bash
#SBATCH -n 32
#SBATCH --ntasks-per-node=32
#SBATCH -p cpu
#SBATCH --mem-per-cpu=3000mb
#SBATCH --exclude work0,work1

source /etc/profile.d/modules.sh

module list
module load vasp/6.2.0-aocl 
module list

ulimit -s unlimited

mpirun --version
echo "#############"
echo "##LIBRARIES##"
echo "#############"
echo $LD_LIBRARY_PATH
echo "#############"
echo "#EXECUTABLE##"
echo "#############"
which vasp_std
echo "#############"
echo "###WRK DIR###"
echo "#############"
pwd

hostname 
export OMP_NUM_THREADS=1

mpirun --bind-to none vasp_std > vasp.out


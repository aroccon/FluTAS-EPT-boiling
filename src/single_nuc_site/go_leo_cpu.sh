#!/bin/bash
#SBATCH --account="IscrB_BMO"
#SBATCH --job-name="boiling_cpu"
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --output=test_cpu.out
#SBATCH -p boost_usr_prod
#SBATCH --error=test_cpu.err
#SBATCH --qos=boost_qos_dbg

# CPU reference run of the same case (binary built with ../compile_cpu.sh)
# same modules used in ../compile_cpu.sh
module purge
module load profile/candidate
module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0-cuda-12.2
export LD_LIBRARY_PATH=/leonardo/prod/spack/06/install/0.22/linux-rhel8-icelake/gcc-12.2.0/fftw-3.3.10-f3ms2lmqzjugb43xzrxjkvmvohumyruy/lib:$LD_LIBRARY_PATH

# debug: print NaN/Inf counts of the main fields after each stage of the time step
#export FLUTAS_NANCHECK=1

SRC_DIR="$(dirname "$(pwd)")"
echo "Using executable: $SRC_DIR/flutas"
# with -D_DECOMP_Z (GNU target) dims_in = 1 4 splits y in 4 (80/4 = 20)
mpirun -np 4 --map-by node:PE=8 $SRC_DIR/flutas

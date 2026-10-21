#!/bin/bash
#SBATCH --account="IscrB_BMO"
#SBATCH --job-name="boiling"
#SBATCH --time=00:30:00
#SBATCH --nodes=1      ##adjust
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=8
#SBATCH --output=test.out
#SBATCH -p boost_usr_prod
#SBATCH --error=test.err
#SBATCH --qos=boost_qos_dbg

# same modules used in ../compile_leo_gpu.sh
module load profile/candidate
module load nvhpc/25.3
module load hpcx-mpi/2.19
# no fftw module: it brings in another MPI (mpirun from nvhpc-24.5) that clashes with hpcx-mpi/2.19;
# libfftw3 is found through the rpath set in targets/target.generic-gpu
echo "Using mpirun: $(which mpirun)"

# executable and GPU/NIC binder are in the src folder (one level up)
SRC_DIR="$(dirname "$(pwd)")"
echo "Using executable: $SRC_DIR/flutas"

# debug: print NaN/Inf counts of the main fields after each stage of the time step
# (see chk_fields in debug.f90); comment out for production runs
#export FLUTAS_NANCHECK=1

# number of MPI ranks (1 per GPU) must match dims_in(1)*dims_in(2) in dns.in
mpirun -np 4 --map-by node:PE=8 --rank-by core $SRC_DIR/binder_leo.sh $SRC_DIR/flutas

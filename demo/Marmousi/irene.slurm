#!/bin/bash
#MSUB -r Marmousi
#MSUB -o Marmousi_%j
#MSUB -e Marmousi_%j.err
#MSUB -N 3
#MSUB -c 10
#-#MSUB -n 128
#MSUB -x
#MSUB -T 10000
#-#MSUB -T 1000
#MSUB -A gen10596
#MSUB -q rome
#MSUB -m store,work,scratch

#module purge
#module load intel/20.0.0 mpi/openmpi

source $ALL_CCCWORKDIR/fftw3.mumps.5.4.0.intel.20.rh8.module

ulimit -s unlimited
ulimit -c unlimited

echo
echo SLURM_NTASKS        $SLURM_NTASKS
echo SLURM_JOB_NUM_NODES $SLURM_JOB_NUM_NODES
echo SLURM_CPUS_PER_TASK $SLURM_CPUS_PER_TASK
echo

# number of threads OpenMP
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# binding OpenMP
export OMP_PLACES=cores

# generate observed data at original datum

cp redatum_obs_org.par redatum.par
time mpirun -n  $SLURM_NTASKS ../../bin/main_redatum_DAW > out1

mv dataobs_00001.bin dataobs_org_00001.bin

# compute DAWs

cp redatum_DAW.par redatum.par
time mpirun -n  $SLURM_NTASKS ../../bin/main_redatum_DAW > out2

# resample the DAW at new datum

cp redatum_Resample.par redatum.par
time mpirun -n  $SLURM_NTASKS ../../bin/main_redatum_DAW > out3

# generate observed data at new datum

cp redatum_obs_new.par redatum.par
time mpirun -n  $SLURM_NTASKS ../../bin/main_redatum_DAW > out4

mv dataobs_00001.bin dataobs_new_00001.bin

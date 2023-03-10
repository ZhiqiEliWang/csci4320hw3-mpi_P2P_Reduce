#!/bin/bash -x

module load spectrum-mpi

#####################################################################################################
# Launch N tasks per compute node allocated. Per below this launches 32 MPI rank per compute node.
# taskset insures that hyperthreaded cores are skipped.
#####################################################################################################
# taskset -c 0-159:4 mpirun -N 32 /gpfs/u/home/SPNR/SPNRcaro/scratch/MPI-Examples/mpi-hello
#mpicc` comp_reduction.c -o comp_reduction`
taskset -c 0-159:4 mpirun -N 32 /gpfs/u/home/PCPC/PCPCwnww/scratch/csci4320hw3/comp_reduction

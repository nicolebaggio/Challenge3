#!/bin/bash

# Number of processors to test
PROCS=(1 2 4)

# Grid size
N=64

# Output data folder
DATA_FOLDER="data"
mkdir -p $DATA_FOLDER  # Create the data folder if it doesn't exist

# Compile the program
mpicxx -fopenmp -o hybrid_parallel_program main.cpp

# Run tests with different number of processors
for PROC in "${PROCS[@]}"
do
  export OMP_NUM_THREADS=$PROC
  mpirun -np $PROC ./hybrid_parallel_program << EOF
$N
$PROC
EOF
done


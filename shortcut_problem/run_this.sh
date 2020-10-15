#! /bin/bash

echo "Removing previous files ..."
rm -f *.o *.mod *.s main

echo "Building the executable ..."

# Add -S for generating assembly code
gfortran -O3 -march=native -fopenmp -c mod_step.f90
gfortran -O3 -march=native -fopenmp main.f90 mod_step.o -o main

echo "Running the executable ..."
export OMP_NUM_THREADS=2
./main
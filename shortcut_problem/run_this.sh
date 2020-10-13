#! /bin/bash

echo "Removing previous files ..."
rm -f *.o *.mod *.s main

echo "Building the executable ..."

# Add -S for generating assembly code
gfortran -c -O3 -march=native mod_step.f90
gfortran -c -O3 -march=native main.f90
gfortran -O3 -march=native *.o -o main

echo "Running the executable ..."
./main
#!/bin/bash

module load intel/2024a
module load CMake/3.29.3-GCCcore-13.3.0

. env/dependencies/install.suitesparse.sh intel icx ifx
. env/dependencies/install.gklib.sh intel icx
. env/dependencies/install.metis32.sh intel icx
. env/dependencies/install.parmetis32.sh intel mpiicx

export CXX=mpiicpx
export CXXFLAGS+=" -fmax-errors=1"

export ESPRESO_USE_WRAPPER_DNBLAS=mkl
export ESPRESO_USE_WRAPPER_DNSOLVER=mkl
export ESPRESO_USE_WRAPPER_LAPACK=mkl
export ESPRESO_USE_WRAPPER_SPBLAS=mkl
export ESPRESO_USE_WRAPPER_SPSOLVER=suitesparse
export ESPRESO_USE_WRAPPER_SCSOLVER=scsolvertriangular
export ESPRESO_USE_WRAPPER_GPU=empty



# export OMP_NUM_THREADS=XXX
# mpirun -n YYY ./build/espreso -c espreso.ecf

#!/bin/bash



module load LUMI/24.03
module load partition/G

if [ ! -d "dependencies/oneAPI-2024.2.0/install" ]
then
    echo "intel toolkit not installed in dependencies, please install it manually"
    return 1
fi
source dependencies/oneAPI-2024.2.0/install/setvars.sh



. dependencies/install.cmake.sh x86_64
. dependencies/install.gklib.sh icx
. dependencies/install.metis32.sh icx
. dependencies/install.parmetis32.sh mpiicx



export CXX=mpiicpx
export CXXFLAGS+=" -ferror-limit=1"

export OMP_PLACES=cores
export OMP_PROC_BIND=close
# export OMP_NUM_THREADS=7,1

export ESPRESO_USE_WRAPPER_DNBLAS=mkl
export ESPRESO_USE_WRAPPER_DNSOLVER=mkl
export ESPRESO_USE_WRAPPER_LAPACK=mkl
export ESPRESO_USE_WRAPPER_SPBLAS=mkl
export ESPRESO_USE_WRAPPER_SPSOLVER=mkl
export ESPRESO_USE_WRAPPER_SCSOLVER=mkl

#!/bin/bash



module load LUMI/24.03
module load partition/G

if [ ! -d "dependencies/oneAPI-2024.2.0/install" ]
then
    echo "intel toolkit not installed in dependencies, please install it manually"
    return 1
fi
source dependencies/oneAPI-2024.2.0/install/setvars.sh



. env/dependencies/install.cmake.sh intel x86_64
. env/dependencies/install.gklib.sh intel icx
. env/dependencies/install.metis32.sh intel icx
. env/dependencies/install.parmetis32.sh intel mpiicx



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

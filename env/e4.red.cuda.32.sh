#!/bin/bash

if [[ "$(hostname)" != *"ngnode"* ]]
then
    echo "ERROR: You have to be on the grace-hopper compute node."
    return
fi



ml purge
ml nvidia/nvhpc-23.11

export CPATH="/opt/share/libs/nvidia/hpc_sdk//Linux_aarch64/23.11/cuda/include:${CPATH}"
export LIBRARY_PATH="/opt/share/libs/nvidia/hpc_sdk//Linux_aarch64/23.11/cuda/lib64:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="/opt/share/libs/nvidia/hpc_sdk//Linux_aarch64/23.11/cuda/lib64:${LD_LIBRARY_PATH}"

export CPATH="/opt/share/libs/nvidia/hpc_sdk/Linux_aarch64/23.11/compilers/include/lp64:${CPATH}"
export LIBRARY_PATH="/opt/share/libs/nvidia/hpc_sdk/Linux_aarch64/23.11/compilers/lib:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="/opt/share/libs/nvidia/hpc_sdk/Linux_aarch64/23.11/compilers/lib:${LD_LIBRARY_PATH}"

export LIBRARY_PATH="/opt/share/libs/nvidia/hpc_sdk//Linux_aarch64/23.11/math_libs/lib64:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="/opt/share/libs/nvidia/hpc_sdk//Linux_aarch64/23.11/math_libs/lib64:${LD_LIBRARY_PATH}"



. env/dependencies/install.cmake.sh aarch64
. env/dependencies/install.suitesparse.sh cuda gcc gfortran "-DBLAS_LIBRARIES=libblas.so -DLAPACK_LIBRARIES=liblapack.so"
. env/dependencies/install.gklib.sh cuda gcc
. env/dependencies/install.metis32.sh cuda gcc
. env/dependencies/install.parmetis32.sh cuda mpicc



export CXX=mpic++
export CXXFLAGS+=" -fmax-errors=1"

export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=12,1

export ESPRESO_RANK_TO_GPU_MAP="0"
export ESPRESO_FORBID_CHECK_SIMD="1" # temporary workaround

export ESPRESO_USE_WRAPPER_DNBLAS=blas
export ESPRESO_USE_WRAPPER_DNSOLVER=lapack
export ESPRESO_USE_WRAPPER_LAPACK=lapack
export ESPRESO_USE_WRAPPER_SPBLAS=suitesparse
export ESPRESO_USE_WRAPPER_SPSOLVER=suitesparse
export ESPRESO_USE_WRAPPER_SCSOLVER=suitesparse
export ESPRESO_USE_WRAPPER_GPU=cuda



# mpirun -n 1 --bind-to numa ./build/espreso -c espreso.ecf

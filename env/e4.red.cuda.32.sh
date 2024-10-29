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



. dependencies/install.cmake.sh
. dependencies/install.suitesparse.sh gcc gfortran "-DBLAS_LIBRARIES=libblas.so -DLAPACK_LIBRARIES=liblapack.so"
. dependencies/install.gklib.sh gcc
. dependencies/install.metis32.sh gcc
. dependencies/install.parmetis32.sh mpicc



export CXX=mpic++
export CXXFLAGS+=" -fmax-errors=1"

export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=12,1

export ESPRESO_RANK_TO_GPU_MAP="0"
export ESPRESO_FORBID_CHECK_SIMD="1" # temporary workaround



# mpirun -n 1 --bind-to numa ./build/espreso -c espreso.ecf

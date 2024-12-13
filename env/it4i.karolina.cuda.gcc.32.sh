#!/bin/bash

cusparse_version="legacy"
# cusparse_version="modern"
if [ "$#" -ge 1 ]
then
    cusparse_version="$1"
fi



ml LAPACK/3.10.0-GCC-12.2.0
ml GCC/11.3.0
ml OpenBLAS/0.3.20-GCC-11.3.0
ml OpenMPI/4.1.4-GCC-11.3.0
ml CMake/3.24.3-GCCcore-11.3.0
if [ "${cusparse_version}" = "legacy" ]; then ml CUDA/11.7.0; fi
if [ "${cusparse_version}" = "modern" ]; then ml CUDA/12.4.0; fi



. dependencies/install.suitesparse.sh gcc gfortran
. dependencies/install.gklib.sh gcc
. dependencies/install.metis32.sh gcc
. dependencies/install.parmetis32.sh mpicc



export CXX=mpic++
export ES_INT_WIDTH=32
export CXXFLAGS+=" -fmax-errors=1"
export BLAS_LIBRARIES=openblas

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=16,1

export ESPRESO_RANK_TO_GPU_MAP="2,3,0,1,6,7,4,5"

if [ "${cusparse_version}" = "legacy" ]; then
    export ESPRESO_USE_CUSPARSE_LEGACY="true"
fi

export ESPRESO_USE_WRAPPER_DNBLAS=blas
export ESPRESO_USE_WRAPPER_DNSOLVER=lapack
export ESPRESO_USE_WRAPPER_LAPACK=lapack
export ESPRESO_USE_WRAPPER_SPBLAS=suitesparse
export ESPRESO_USE_WRAPPER_SPSOLVER=suitesparse
export ESPRESO_USE_WRAPPER_SCSOLVER=suitesparse
export ESPRESO_USE_WRAPPER_GPU=cuda



# mpirun -n 8 --bind-to numa ./build/espreso -c espreso.ecf

#!/bin/bash



ml CMake/3.29.3-GCCcore-13.3.0
ml intel/2024a

. env/dependencies/install.gklib.sh cudamkl icx
. env/dependencies/install.metis32.sh cudamkl icx
. env/dependencies/install.parmetis32.sh cudamkl mpiicx
. env/dependencies/install.suitesparse.sh cudamkl icx ifx

ml CUDA/12.8.0



export CXX=mpiicpx
export ES_INT_WIDTH=32
export CXXFLAGS+=" -fmax-errors=1"
export CXXFLAGS+=" -DESPRESO_STACKTIMER_ENABLE"

export ESPRESO_CUDA_ALLOW_UNSUPPORTED_COMPILER=1

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=16,1

export ESPRESO_RANK_TO_GPU_MAP="2,3,0,1,6,7,4,5"

export ESPRESO_USE_WRAPPER_DNBLAS=mkl
export ESPRESO_USE_WRAPPER_DNSOLVER=mkl
export ESPRESO_USE_WRAPPER_LAPACK=mkl
export ESPRESO_USE_WRAPPER_SPBLAS=mkl
export ESPRESO_USE_WRAPPER_SPSOLVER=suitesparse
export ESPRESO_USE_WRAPPER_SCSOLVER=mkl
export ESPRESO_USE_WRAPPER_GPU=cuda

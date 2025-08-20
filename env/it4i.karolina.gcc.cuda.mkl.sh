#!/bin/bash

buildname="build-it4i-karolina-gcc-cuda-mkl"
export WAFLOCK=".lock-waf_linux_${buildname}"
curr_buildname="$(readlink build)"
if [ "${curr_buildname}" != "${buildname}" ]
then
    rm -rf build
    ln -s "${buildname}" build
fi



ml CMake/3.29.3-GCCcore-13.3.0
ml intel/2024a
ml GCC/13.3.0
ml CUDA/12.4.0
ml OpenMPI/5.0.3-GCC-13.3.0



. env/dependencies/install.gklib.sh cudamkl gcc
. env/dependencies/install.metis32.sh cudamkl gcc
. env/dependencies/install.parmetis32.sh cudamkl mpicc



export CXX=mpic++
export ES_INT_WIDTH=32
export CXXFLAGS+=" -fmax-errors=1"

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=16,1

export ESPRESO_RANK_TO_GPU_MAP="2,3,0,1,6,7,4,5"

export ESPRESO_USE_WRAPPER_DNBLAS=mkl
export ESPRESO_USE_WRAPPER_DNSOLVER=mkl
export ESPRESO_USE_WRAPPER_LAPACK=mkl
export ESPRESO_USE_WRAPPER_SPBLAS=mkl
export ESPRESO_USE_WRAPPER_SPSOLVER=mkl
export ESPRESO_USE_WRAPPER_SCSOLVER=mkl
export ESPRESO_USE_WRAPPER_GPU=cuda



# mpirun -n 8 --bind-to numa ./build/espreso -c espreso.ecf

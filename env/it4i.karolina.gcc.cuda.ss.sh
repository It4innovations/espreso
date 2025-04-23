#!/bin/bash

cudaversionname="${1}"
cudaversion=""
if [ "${cudaversionname}" == "legacy" ]
then
    cudaversion="11.7.0"
fi
if [ "${cudaversionname}" == "modern" ]
then
    cudaversion="12.8.0"
fi
if [ "${cudaversion}" == "" ]
then
    echo "wrong argument - cuda version name - legacy/modern"
    return
fi



buildname="build-it4i-karolina-gcc-cuda${cudaversionname}-ss"
export WAFLOCK=".lock-waf_linux_${buildname}"
rm -rf build
ln -s "${buildname}" build



ml LAPACK/3.10.0-GCC-12.2.0
ml GCC/11.3.0
ml OpenBLAS/0.3.20-GCC-11.3.0
ml OpenMPI/4.1.4-GCC-11.3.0
ml CMake/3.24.3-GCCcore-11.3.0
ml "CUDA/${cudaversion}"



. env/dependencies/install.suitesparse.sh gcccudass gcc gfortran
. env/dependencies/install.gklib.sh gcccudass gcc
. env/dependencies/install.metis32.sh gcccudass gcc
. env/dependencies/install.parmetis32.sh gcccudass mpicc



export CXX=mpic++
export ES_INT_WIDTH=32
export CXXFLAGS+=" -fmax-errors=1"
export BLAS_LIBRARIES=openblas

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=16,1

export ESPRESO_RANK_TO_GPU_MAP="2,3,0,1,6,7,4,5"

if [ "${cudaversionname}" = "legacy" ]; then
    export ESPRESO_USE_CUSPARSE_LEGACY="true"
fi

export ESPRESO_USE_WRAPPER_DNBLAS=blas
export ESPRESO_USE_WRAPPER_DNSOLVER=lapack
export ESPRESO_USE_WRAPPER_LAPACK=lapack
export ESPRESO_USE_WRAPPER_SPBLAS=suitesparse
export ESPRESO_USE_WRAPPER_SPSOLVER=suitesparse
export ESPRESO_USE_WRAPPER_SCSOLVER=suitesparse
export ESPRESO_USE_WRAPPER_GPU=cuda

export CXXFLAGS+=" -DESPRESO_STACKTIMER_ENABLE"



# mpirun -n 8 --bind-to numa ./build/espreso -c espreso.ecf

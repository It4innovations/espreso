#!/bin/bash

cudaversionname="${1}"
cudaversion=""
if [ "${cudaversionname}" == "legacy" ]
then
    cudaversion="11.8.0"
fi
if [ "${cudaversionname}" == "modern" ]
then
    cudaversion="12.8.0"
fi
if [ "${cudaversion}" == "" ]
then
    echo "mising argument - cuda version name - legacy/modern"
    return
fi

ml gcc/11.4.0
ml cmake/3.30.5
ml openmpi/4.1.5-gcc
ml mkl/2024.2

. env/dependencies/install.gklib.sh gcccudamklss gcc
. env/dependencies/install.metis32.sh gcccudamklss gcc
. env/dependencies/install.parmetis32.sh gcccudamklss mpicc
. env/dependencies/install.suitesparse.sh gcccudamklss gcc gfortran

export CUDA_ROOT="${PWD}/dependencies/cuda-${cudaversion}/install"
if [ ! -d "${CUDA_ROOT}" ]
then
    echo "please download and install cuda-${cudaversion}"
    # just make symlink to the actual install dir, or
    # mkdir -p "dependencies/cuda-12.8.0"
    # cd "dependencies/cuda-12.8.0"
    # wget https://developer.download.nvidia.com/compute/cuda/12.8.0/local_installers/cuda_12.8.0_570.86.10_linux.run
    # sh cuda_12.8.0_570.86.10_linux.run --silent --toolkit --installpath=${PWD}/install
fi
export PATH="${CUDA_ROOT}/bin:${PATH}"
export CPATH="${CUDA_ROOT}/include:${CPATH}"
export LIBRARY_PATH="${CUDA_ROOT}/lib64:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${CUDA_ROOT}/lib64:${LD_LIBRARY_PATH}"



export CXX=mpic++
export ES_INT_WIDTH=32
export CXXFLAGS+=" -fmax-errors=1"
export CXXFLAGS+=" -DESPRESO_STACKTIMER_ENABLE"

if [ "${cudaversionname}" = "legacy" ]; then
    export ESPRESO_USE_CUSPARSE_LEGACY="true"
fi

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=20,1

export ESPRESO_RANK_TO_GPU_MAP="0,1,2,3"

export ESPRESO_USE_WRAPPER_DNBLAS=mkl
export ESPRESO_USE_WRAPPER_DNSOLVER=mkl
export ESPRESO_USE_WRAPPER_LAPACK=mkl
export ESPRESO_USE_WRAPPER_SPBLAS=mkl
export ESPRESO_USE_WRAPPER_SPSOLVER=suitesparse
export ESPRESO_USE_WRAPPER_SCSOLVER=mkl
export ESPRESO_USE_WRAPPER_GPU=cuda



# mpirun -n 1 --bind-to numa ...

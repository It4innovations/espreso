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



buildname="build-it4i-karolina-gcc-cuda${cudaversionname}-mkl-ss"
export WAFLOCK=".lock-waf_linux_${buildname}"
rm -rf build
ln -s "${buildname}" build



ml intel/2024a
if [ "${cudaversionname}" = "legacy" ]; then
    ml GCC/11.3.0
    ml CMake/3.24.3-GCCcore-11.3.0
    ml OpenMPI/4.1.4-GCC-11.3.0
fi
if [ "${cudaversionname}" = "modern" ]; then
    ml CMake/3.29.3-GCCcore-13.3.0
    ml OpenMPI/5.0.3-GCC-13.3.0
fi

. env/dependencies/install.gklib.sh gcccudamklss gcc
. env/dependencies/install.metis32.sh gcccudamklss gcc
. env/dependencies/install.parmetis32.sh gcccudamklss mpicc
. env/dependencies/install.suitesparse.sh gcccudamklss gcc gfortran
. env/dependencies/install.mumps.sh gcccudamklss mpicc mpifort "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core" "-lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core" "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core"

ml "CUDA/${cudaversion}"



export CXX=mpic++
export ES_INT_WIDTH=32
export CXXFLAGS+=" -fmax-errors=1"
export CXXFLAGS+=" -DESPRESO_STACKTIMER_ENABLE"

if [ "${cudaversionname}" = "legacy" ]; then
    export ESPRESO_USE_CUSPARSE_LEGACY="true"
fi

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



# mpirun -n XYZ --bind-to numa ...

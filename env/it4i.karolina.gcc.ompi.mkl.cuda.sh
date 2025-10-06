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



# build-institution-machine-compiler-mpi-blas-specialname
buildname="build-it4i-karolina-gcc-ompi-mkl-cuda${cudaversionname}"
export WAFLOCK=".lock-waf_linux_${buildname}"
curr_buildname="$(readlink build)"
if [ "${curr_buildname}" != "${buildname}" ]
then
    rm -rf build
    ln -s "${buildname}" build
fi



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

ml "CUDA/${cudaversion}"

. env/dependencies/install.gklib.sh gcc_ompi_mkl_cuda gcc
. env/dependencies/install.metis32.sh gcc_ompi_mkl_cuda gcc
. env/dependencies/install.parmetis32.sh gcc_ompi_mkl_cuda mpicc
. env/dependencies/install.suitesparse.sh gcc_ompi_mkl_cuda gcc gfortran
# . env/dependencies/install.mumps.sh gcc_ompi_mkl_cuda mpicc mpifort "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core" "-lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core" "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core"
# . env/dependencies/install.strumpack.sh gcc_ompi_mkl_cuda g++ gcc gfortran "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core" "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core"
. env/dependencies/install.pastix.sh gcc_ompi_mkl_cuda g++ gcc gfortran Intel10_64lp_seq
# . env/dependencies/install.superlu_dist.sh gcc_ompi_mkl_cuda g++ gcc "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core"



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

ngpus="$(nvidia-smi -L | wc -l)"
if [ "${ngpus}" -eq 8 ]; then
    export ESPRESO_MAP_LOCALRANK_TO_GPU="2 3 0 1 6 7 4 5"
elif [ "${ngpus}" -ne 1 ]; then
    echo
    echo "WARNING: only support 1 or 8 GPUs per node"
    echo "  defaulting to using linear mpi-gpu map"
    echo "  set your own ESPRESO_MAP_LOCALRANK_TO_GPU env var to override"
    echo
fi

export ESPRESO_USE_WRAPPER_DNBLAS=mkl
export ESPRESO_USE_WRAPPER_DNSOLVER=mkl
export ESPRESO_USE_WRAPPER_LAPACK=mkl
export ESPRESO_USE_WRAPPER_SPBLAS=mkl
export ESPRESO_USE_WRAPPER_SPSOLVER=suitesparse
export ESPRESO_USE_WRAPPER_SCSOLVER=mkl
export ESPRESO_USE_WRAPPER_GPU=cuda



# mpirun -n XYZ --bind-to numa ...

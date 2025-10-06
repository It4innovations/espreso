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
    echo "wrong argument - cuda version name - legacy/modern"
    return
fi



# build-institution-machine-compiler-mpi-blas-specialname
buildname="build-bsc-mn5-gcc-ompi-mkl-cuda${cudaversionname}"
export WAFLOCK=".lock-waf_linux_${buildname}"
curr_buildname="$(readlink build)"
if [ "${curr_buildname}" != "${buildname}" ]
then
    rm -rf build
    ln -s "${buildname}" build
fi



ml gcc/11.4.0
ml cmake/3.30.5
ml openmpi/4.1.5-gcc
ml mkl/2024.2

. env/dependencies/install.gklib.sh gcc_ompi_mkl_cuda gcc
. env/dependencies/install.metis32.sh gcc_ompi_mkl_cuda gcc
. env/dependencies/install.parmetis32.sh gcc_ompi_mkl_cuda mpicc
. env/dependencies/install.suitesparse.sh gcc_ompi_mkl_cuda gcc gfortran
. env/dependencies/install.pastix.sh gcc_ompi_mkl_cuda g++ gcc gfortran Intel10_64lp_seq

export CUDA_ROOT="${PWD}/dependencies/cuda-${cudaversion}"
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

export SRUN_CPUS_PER_TASK=20

# export ESPRESO_MAP_LOCALRANK_TO_GPU="0"

export ESPRESO_USE_WRAPPER_DNBLAS=mkl
export ESPRESO_USE_WRAPPER_DNSOLVER=mkl
export ESPRESO_USE_WRAPPER_LAPACK=mkl
export ESPRESO_USE_WRAPPER_SPBLAS=mkl
export ESPRESO_USE_WRAPPER_SPSOLVER=suitesparse
export ESPRESO_USE_WRAPPER_SCSOLVER=mkl
export ESPRESO_USE_WRAPPER_GPU=cuda



# salloc -A ... -q acc_ehpc -N 1 --ntasks-per-node 4 -c 20 --gres=gpu:4 --threads-per-core 1 -t 1:00:00
# srun -n 1 -c 20 --gpus-per-task 1 ...

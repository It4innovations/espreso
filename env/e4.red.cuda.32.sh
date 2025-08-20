#!/bin/bash

if [[ "$(hostname)" != *"ngnode"* ]]
then
    echo "ERROR: You have to be on the grace-hopper compute node."
    return 1
fi

if [ "$#" -ge 1 ]
then
    cusparse_version="$1"
else
    echo "You have to specify cusparse version (legacy/modern)"
    return 2
fi



buildname="build-e4-red-cuda${cusparse_version}-32"
export WAFLOCK=".lock-waf_linux_${buildname}"
curr_buildname="$(readlink build)"
if [ "${curr_buildname}" != "${buildname}" ]
then
    rm -rf build
    ln -s "${buildname}" build
fi



ml purge

if [ "${cusparse_version}" == "legacy" ]; then cuda_version="11.7.1"; fi
if [ "${cusparse_version}" == "modern" ]; then cuda_version="12.5.1"; fi

cuda_root="${PWD}/dependencies/cuda-${cuda_version}"
if [ ! -d "${cuda_root}" ]
then
    echo "CUDA ${cuda_version} not installed in dependencies, please install it manually"
    # for example:
    # mkdir -p dependencies/cuda-11.7.1
    # cd dependencies/cuda-11.7.1
    # wget https://developer.download.nvidia.com/compute/cuda/11.7.1/local_installers/cuda_11.7.1_515.65.01_linux_sbsa.run
    # sh cuda_11.7.1_515.65.01_linux_sbsa.run --silent --toolkit --installpath=${PWD}
    return 3
fi

# for CUDA compiler and libs
export PATH="${cuda_root}/bin:${PATH}"
export LIBRARY_PATH="${cuda_root}/lib64:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${cuda_root}/lib64:${LD_LIBRARY_PATH}"
export CPATH="${cuda_root}/include:${CPATH}"

openmpi_root="${PWD}/dependencies/openmpi-5.0.6/install"
if [ ! -d "${cuda_root}" ]
then
    echo "OpenMPI 5.0.6 not installed in dependencies, please install it manually"
    # for example:
    # cd dependencies
    # wget https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.6.tar.gz
    # cd openmpi-5.0.6
    # ./configure --prefix="${PWD}/install"
    # make -j$(nproc)
    # make install
    return 3
fi

# openmpi
export PATH="${openmpi_root}/bin:${PATH}"
export LIBRARY_PATH="${openmpi_root}/lib:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${openmpi_root}/lib:${LD_LIBRARY_PATH}"
export CPATH="${openmpi_root}/include:${CPATH}"

# blas and lapack
export CPATH="/opt/share/libs/nvidia/hpc_sdk/Linux_aarch64/23.11/compilers/include/lp64:${CPATH}"
export LIBRARY_PATH="/opt/share/libs/nvidia/hpc_sdk/Linux_aarch64/23.11/compilers/lib:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="/opt/share/libs/nvidia/hpc_sdk/Linux_aarch64/23.11/compilers/lib:${LD_LIBRARY_PATH}"




. env/dependencies/install.cmake.sh cuda aarch64
. env/dependencies/install.suitesparse.sh cuda gcc gfortran "-DBLAS_LIBRARIES=libblas.so -DLAPACK_LIBRARIES=liblapack.so"
. env/dependencies/install.gklib.sh cuda gcc
. env/dependencies/install.metis32.sh cuda gcc
. env/dependencies/install.parmetis32.sh cuda mpicc



if [ "${cusparse_version}" == "legacy" ]
then
    # 90 not supported with 11.7 cuda toolkit
    export CUDAARCH=compute_80
fi

export CXX=mpic++
export CXXFLAGS+=" -fmax-errors=1"
export LINKFLAGS+=" -Wl,--no-as-needed -lnvf"

export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=12,1

if [ "${cusparse_version}" == "legacy" ];
then
    export ESPRESO_USE_CUSPARSE_LEGACY="1"
fi
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

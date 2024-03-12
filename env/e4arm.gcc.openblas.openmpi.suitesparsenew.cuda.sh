#!/bin/bash

if [[ "$(hostname)" != *"ngnode"* ]]
then
    echo "ERROR: You have to be on the grace-hopper compute node."
    return
fi

# cusparse_version="legacy" # legacy cusparse not available, because cuda<12 not available
cusparse_version="modern"



DEPENDENCIES_DIR="dependencies"
if [ ! -d "${DEPENDENCIES_DIR}" ]; then mkdir "${DEPENDENCIES_DIR}"; fi

export MY_CMAKE_VERSION="3.28.3-linux-aarch64"
export MY_CMAKE_DIR="cmake-${MY_CMAKE_VERSION}"
export MY_CMAKE_ROOT="${DEPENDENCIES_DIR}/${MY_CMAKE_DIR}"
if [ ! -d "${MY_CMAKE_ROOT}" ]
then
    (
        echo "Downloading newer CMake"
        cd "${DEPENDENCIES_DIR}"
        wget https://github.com/Kitware/CMake/releases/download/v3.28.3/cmake-3.28.3-linux-aarch64.tar.gz
        tar -xf cmake-3.28.3-linux-aarch64.tar.gz
        rm cmake-3.28.3-linux-aarch64.tar.gz
    )
fi
export PATH="${PWD}/${MY_CMAKE_ROOT}/bin:${PATH}"

# (open)blas is already in /usr/lib64
# lapacke is already in /usr/lib64

export MY_VERSION_METIS="5.1.0"
export MY_METIS_DIR="metis-${MY_VERSION_METIS}"
export MY_METIS_ROOT="${DEPENDENCIES_DIR}/${MY_METIS_DIR}"
if [ ! -d "${MY_METIS_ROOT}" ]
then
    (
        echo "Metis not found, installing locally ..."
        cd "${DEPENDENCIES_DIR}"
        wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
        tar -xf metis-5.1.0.tar.gz
        rm metis-5.1.0.tar.gz
        cd metis-5.1.0
        ml gcc-8.5.0/ompi-4.1.4
        sed -i 's/#define REALTYPEWIDTH 32/#define REALTYPEWIDTH 64/' include/metis.h
        make config shared=1 cc=gcc gklib_path="${PWD}/GKlib" prefix="${PWD}"
        make -j$(nproc)
        make install
    )
fi

VERSION_SUITESPARSE=v7.6.0
export SUITESPARSE_DIR="SuiteSparse_${VERSION_SUITESPARSE}_gcc_openblas"
export SUITESPARSE_ROOT="${DEPENDENCIES_DIR}/${SUITESPARSE_DIR}"
if [ ! -d "${SUITESPARSE_ROOT}" ]
then
    (
        echo "SuiteSparse ${VERSION_SUITESPARSE} not installed, installing locally..."
        cd "${DEPENDENCIES_DIR}"
        git clone -b "${VERSION_SUITESPARSE}" https://github.com/DrTimothyAldenDavis/SuiteSparse.git "${SUITESPARSE_DIR}"
        cd "${SUITESPARSE_DIR}"

        ml gcc-8.5.0/ompi-4.1.4
        mkdir -p build
        cd build
        cmake -DCMAKE_C_COMPILER=gcc -DENABLE_CUDA=false -DSUITESPARSE_ENABLE_PROJECTS="cholmod" -DCMAKE_INSTALL_PREFIX=.. ..
        cmake --build . -j $(nproc)
        cmake --install .
    )
fi



ml gcc-8.5.0/ompi-4.1.4

ml nvidia/cuda-12.3



export CPATH="${PWD}/${SUITESPARSE_ROOT}/include:${CPATH}"
export LIBRARY_PATH="${PWD}/${SUITESPARSE_ROOT}/lib:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${PWD}/${SUITESPARSE_ROOT}/lib:${LD_LIBRARY_PATH}"

export CPATH="${PWD}/${MY_METIS_ROOT}/include:${CPATH}"
export LIBRARY_PATH="${PWD}/${MY_METIS_ROOT}/lib:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH="${PWD}/${MY_METIS_ROOT}/lib:${LD_LIBRARY_PATH}"

export CPATH="/opt/share/libs/nvidia/cuda-12.3/include:${CPATH}"
export LIBRARY_PATH="/opt/share/libs/nvidia/cuda-12.3/lib64:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="/opt/share/libs/nvidia/cuda-12.3/lib64:${LD_LIBRARY_PATH}"

export CXX=mpic++
export CXXFLAGS+=" -fmax-errors=1"

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=72,1

export ESPRESO_RANK_TO_GPU_MAP="0"

if [ "${cusparse_version}" = "legacy" ]; then
    export ESPRESO_USE_CUSPARSE_LEGACY="true"
fi



# mpirun -n 1 ./build/espreso -c espreso.ecf

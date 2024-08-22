#!/bin/bash

DEPENDENCIES_DIR="dependencies"
mkdir -p "${DEPENDENCIES_DIR}"



ml intel/2024a



METIS_DIR="metis-5.1.0_intel2023b"
METIS_ROOT="${DEPENDENCIES_DIR}/${METIS_DIR}"
if [ ! -d "${METIS_ROOT}" ]
then
    (
        ml CMake/3.27.6-GCCcore-13.2.0

        echo "METIS not installed, installing locally..."
        cd "${DEPENDENCIES_DIR}"
        mkdir "${METIS_DIR}"
        cd "${METIS_DIR}"
        # wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
        wget https://github.com/xijunke/METIS-1/raw/master/metis-5.1.0.tar.gz
        tar -xf metis-5.1.0.tar.gz
        rm metis-5.1.0.tar.gz
        mv metis-5.1.0/* .
        rm -r metis-5.1.0
        sed -i 's/#define REALTYPEWIDTH 32/#define REALTYPEWIDTH 64/' include/metis.h
        make config shared=1 cc=icx gklib_path="${PWD}/GKlib" prefix="${PWD}"
        make -j$(nproc)
        make install
    )
fi
export CPATH="${PWD}/${METIS_ROOT}/include:${CPATH}"
export LIBRARY_PATH="${PWD}/${METIS_ROOT}/lib:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH="${PWD}/${METIS_ROOT}/lib:${LD_LIBRARY_PATH}"



VERSION_SUITESPARSE=v7.6.0
export SUITESPARSE_DIR="SuiteSparse_${VERSION_SUITESPARSE}_intel2023b_debug"
export SUITESPARSE_ROOT="${DEPENDENCIES_DIR}/${SUITESPARSE_DIR}"
if [ ! -d "${SUITESPARSE_ROOT}" ]
then
    (
        ml CMake/3.27.6-GCCcore-13.2.0

        echo "SuiteSparse ${VERSION_SUITESPARSE} not installed, installing locally..."
        cd "${DEPENDENCIES_DIR}"
        git clone -b "${VERSION_SUITESPARSE}" https://github.com/DrTimothyAldenDavis/SuiteSparse.git "${SUITESPARSE_DIR}"
        cd "${SUITESPARSE_DIR}"

        mkdir -p build
        cd build
        cmake -DCMAKE_C_COMPILER=icx -DCMAKE_Fortran_COMPILER=ifx -DENABLE_CUDA=false -DSUITESPARSE_ENABLE_PROJECTS="cholmod" -DCMAKE_INSTALL_PREFIX=.. ..
        # cmake -DCMAKE_C_COMPILER=icx -DCMAKE_Fortran_COMPILER=ifx -DSUITESPARSE_USE_64BIT_BLAS=true -DENABLE_CUDA=false -DSUITESPARSE_ENABLE_PROJECTS="cholmod" -DCMAKE_INSTALL_PREFIX=.. ..
        cmake --build . -j $(nproc)
        cmake --install .
    )
fi
export CPATH="${PWD}/${SUITESPARSE_ROOT}/include:${CPATH}"
export LIBRARY_PATH="${PWD}/${SUITESPARSE_ROOT}/lib64:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${PWD}/${SUITESPARSE_ROOT}/lib64:${LD_LIBRARY_PATH}"



export CXX=mpiicpx
export CXXFLAGS+=" -ferror-limit=1"



export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=16,1

export ESPRESO_FORBID_MKL_PARDISO="1"
export ESPRESO_RANK_TO_GPU_MAP="0"

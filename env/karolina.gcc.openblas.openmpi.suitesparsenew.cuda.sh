#!/bin/bash

cusparse_version="legacy"
# cusparse_version="modern"
if [ "$#" -ge 1 ]
then
    cusparse_version="$1"
fi



DEPENDENCIES_DIR="dependencies"
if [ ! -d "${DEPENDENCIES_DIR}" ]; then mkdir "${DEPENDENCIES_DIR}"; fi

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

        ml LAPACK/3.10.0-GCC-12.2.0
        ml GCC/11.3.0
        ml OpenBLAS/0.3.20-GCC-11.3.0
        ml CMake/3.24.3-GCCcore-11.3.0
        mkdir -p build
        cd build
        cmake -DCMAKE_C_COMPILER=gcc -DENABLE_CUDA=false -DSUITESPARSE_ENABLE_PROJECTS="cholmod" -DCMAKE_INSTALL_PREFIX=.. ..
        cmake --build . -j $(nproc)
        cmake --install .
    )
fi
export CPATH="${PWD}/${SUITESPARSE_ROOT}/include:${CPATH}"
export LIBRARY_PATH="${PWD}/${SUITESPARSE_ROOT}/lib64:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${PWD}/${SUITESPARSE_ROOT}/lib64:${LD_LIBRARY_PATH}"



ml LAPACK/3.10.0-GCC-12.2.0
ml GCC/11.3.0
ml OpenBLAS/0.3.20-GCC-11.3.0
ml OpenMPI/4.1.4-GCC-11.3.0
ml METIS/5.1.0-GCCcore-11.3.0
if [ "${cusparse_version}" = "legacy" ]; then ml CUDA/11.7.0; fi
if [ "${cusparse_version}" = "modern" ]; then ml CUDA/12.4.0; fi

export CXX=mpic++
export ES_INT_WIDTH=32
export CXXFLAGS+=" -fmax-errors=1"
export BLAS_LIBRARIES=openblas

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=16,1

export ESPRESO_RANK_TO_GPU_MAP="2,3,0,1,6,7,4,5"

if [ "${cusparse_version}" = "legacy" ]; then
    export ESPRESO_USE_CUSPARSE_LEGACY="true"
fi



# mpirun -n 8 --bind-to numa ./build/espreso -c espreso.ecf

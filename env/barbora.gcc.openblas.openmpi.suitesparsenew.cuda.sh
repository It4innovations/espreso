#!/bin/bash

cusparse_version="legacy"
# cusparse_version="modern"



DEPENDENCIES_DIR="dependencies"
if [ ! -d "${DEPENDENCIES_DIR}" ]; then mkdir "${DEPENDENCIES_DIR}"; fi

LAPACK_ROOT="${DEPENDENCIES_DIR}/lapack"
if [ ! -d "${LAPACK_ROOT}" ]
then
    (
        ml GCC/11.3.0
        ml OpenBLAS/0.3.20-GCC-11.3.0
        ml CMake/3.24.3-GCCcore-11.3.0
        cd "${DEPENDENCIES_DIR}"
        wget http://www.netlib.org/lapack/lapack.tgz
        tar -xf lapack.tgz
        rm lapack.tgz
        mv lapack-3.10.1 lapack
        cd lapack
        mkdir build
        cd build
        cmake -DCMAKE_C_COMPILER=gcc -DLAPACKE=true -DBLAS_LIBRARIES=/apps/all/OpenBLAS/0.3.20-GCC-11.3.0/lib/libopenblas.so -DCMAKE_INSTALL_PREFIX=.. ..
        make -j $(nproc)
        make install
    )
fi
export CPATH="${PWD}/${LAPACK_ROOT}/include:${CPATH}"
export LIBRARY_PATH="${PWD}/${LAPACK_ROOT}/lib64:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${PWD}/${LAPACK_ROOT}/lib64:${LD_LIBRARY_PATH}"

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



ml METIS/5.1.0-GCCcore-10.3.0
ml GCC/11.3.0
ml OpenBLAS/0.3.20-GCC-11.3.0
ml OpenMPI/4.1.4-GCC-11.3.0
if [ "${cusparse_version}" = "legacy" ]; then ml CUDA/11.7.0; fi
if [ "${cusparse_version}" = "modern" ]; then ml CUDA/12.2.0; fi

export CXX=mpic++
export CXXFLAGS+=" -fmax-errors=1"
export BLAS_LIBRARIES=openblas

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=6,1

export ESPRESO_RANK_TO_GPU_MAP="0,1,2,3"

if [ "${cusparse_version}" = "legacy" ]; then
    export ESPRESO_USE_CUSPARSE_LEGACY="true"
fi



# mpirun -n 4 --map-by node:PE=6 ./build/espreso -c espreso.ecf

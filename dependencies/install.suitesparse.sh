#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"

VERSION_SUITESPARSE=v7.6.0

SUITESPARSE_DIR="SuiteSparse_${VERSION_SUITESPARSE}"
SUITESPARSE_ROOT="${DEPENDENCIES_DIR}/${SUITESPARSE_DIR}"
if [ ! -d "${SUITESPARSE_ROOT}" ]
then
  sh ${DEPENDENCIES_DIR}/clone.suitesparse.sh
fi

if [ ! -d "${SUITESPARSE_ROOT}/$1" ]
then
    (
        cd "${SUITESPARSE_ROOT}"
        mkdir -p build
        cd build
        cmake -DCMAKE_C_COMPILER=$1 -DENABLE_CUDA=false -DSUITESPARSE_ENABLE_PROJECTS="cholmod" -DCMAKE_INSTALL_PREFIX=${SUITESPARSE_ROOT}/$1 $2 ..
        cmake --build . -j $(nproc)
        cmake --install .
    )
fi

export CPATH="${SUITESPARSE_ROOT}/$1/include:${CPATH}"
export LIBRARY_PATH="${SUITESPARSE_ROOT}/$1/lib:${SUITESPARSE_ROOT}/$1/lib64:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${SUITESPARSE_ROOT}/$1/lib:${SUITESPARSE_ROOT}/$1/lib64:${LD_LIBRARY_PATH}"

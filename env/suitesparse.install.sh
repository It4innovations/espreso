#!/bin/bash

DEPENDENCIES_DIR="dependencies"
if [ ! -d "${DEPENDENCIES_DIR}" ]; then mkdir "${DEPENDENCIES_DIR}"; fi

VERSION_SUITESPARSE=v7.6.0

export SUITESPARSE_DIR="SuiteSparse_${VERSION_SUITESPARSE}"
export SUITESPARSE_ROOT="${DEPENDENCIES_DIR}/${SUITESPARSE_DIR}"
if [ ! -d "${SUITESPARSE_ROOT}/build$1" ]
then
    (
        cd "${SUITESPARSE_ROOT}"
        mkdir -p build$1
        cd build$1
        cmake -DCMAKE_C_COMPILER=gcc -DENABLE_CUDA=false -DSUITESPARSE_ENABLE_PROJECTS="cholmod" -DCMAKE_INSTALL_PREFIX=. ..
        cmake --build . -j $(nproc)
        cmake --install .
    )
fi

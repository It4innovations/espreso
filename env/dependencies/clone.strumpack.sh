#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

VERSION_STRUMPACK=v8.0.0

STRUMPACK_DIR="Strumpack_${VERSION_STRUMPACK}"
STRUMPACK_ROOT="${DEPENDENCIES_DIR}/${STRUMPACK_DIR}"
if [ ! -d "${STRUMPACK_ROOT}" ]
then
    (
        echo "Strumpack ${VERSION_STRUMPACK} not found, cloning..."
        cd "${DEPENDENCIES_DIR}"
        git clone -b "${VERSION_STRUMPACK}" https://github.com/pghysels/STRUMPACK.git "${STRUMPACK_DIR}"
        sed -i 's/DenseMW_t X(N, nrhs, x, N);/DenseMW_t X(N, nrhs, x, ldx);/g' "${STRUMPACK_DIR}/src/SparseSolverBase.cpp"
    )
fi

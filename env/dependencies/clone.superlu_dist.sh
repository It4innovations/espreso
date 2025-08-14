#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

VERSION_SUPERLU=v9.1.0

SUPERLU_DIR="SuperLU_dist_${VERSION_SUPERLU}"
SUPERLU_ROOT="${DEPENDENCIES_DIR}/${SUPERLU_DIR}"
if [ ! -d "${SUPERLU_ROOT}" ]
then
    (
        echo "SuperLU_dist ${VERSION_SUPERLU} not found, cloning..."
        cd "${DEPENDENCIES_DIR}"
        git clone -b "${VERSION_SUPERLU}" https://github.com/xiaoyeli/superlu_dist.git "${SUPERLU_DIR}"
        cd "${SUPERLU_DIR}"
        sed -i 's|#pragma omp|// #pragma omp|g' SRC/double/dSchCompUdt-gpu.c
    )
fi

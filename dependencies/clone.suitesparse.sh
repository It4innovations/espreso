#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"

VERSION_SUITESPARSE=v7.6.0

SUITESPARSE_DIR="SuiteSparse_${VERSION_SUITESPARSE}"
SUITESPARSE_ROOT="${DEPENDENCIES_DIR}/${SUITESPARSE_DIR}"
if [ ! -d "${SUITESPARSE_ROOT}" ]
then
    (
        echo "SuiteSparse ${VERSION_SUITESPARSE} not found, cloning..."
        cd "${DEPENDENCIES_DIR}"
        git clone -b "${VERSION_SUITESPARSE}" https://github.com/DrTimothyAldenDavis/SuiteSparse.git "${SUITESPARSE_DIR}"
    )
fi

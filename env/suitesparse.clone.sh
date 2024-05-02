#!/bin/bash

DEPENDENCIES_DIR="dependencies"
if [ ! -d "${DEPENDENCIES_DIR}" ]; then mkdir "${DEPENDENCIES_DIR}"; fi

VERSION_SUITESPARSE=v7.6.0

export SUITESPARSE_DIR="SuiteSparse_${VERSION_SUITESPARSE}"
export SUITESPARSE_ROOT="${DEPENDENCIES_DIR}/${SUITESPARSE_DIR}"
if [ ! -d "${SUITESPARSE_ROOT}" ]
then
    (
        echo "SuiteSparse ${VERSION_SUITESPARSE} not found, cloning..."
        cd "${DEPENDENCIES_DIR}"
        git clone -b "${VERSION_SUITESPARSE}" https://github.com/DrTimothyAldenDavis/SuiteSparse.git "${SUITESPARSE_DIR}"
    )
fi

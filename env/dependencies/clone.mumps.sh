#!/bin/bash

ORIG_DIR="${PWD}"

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

VERSION_MUMPS="${1}"

MUMPS_DIR="MUMPS_${VERSION_MUMPS}"
MUMPS_ROOT="${DEPENDENCIES_DIR}/${MUMPS_DIR}"
if [ ! -d "${MUMPS_ROOT}" ]
then
    (
        echo "MUMPS ${VERSION_MUMPS} not found, cloning..."
        cd "${DEPENDENCIES_DIR}"
        "${ORIG_DIR}/env/dependencies/download.sh" "https://mumps-solver.org/MUMPS_${VERSION_MUMPS}.tar.gz"
        tar -xf "MUMPS_${VERSION_MUMPS}.tar.gz"
        rm "MUMPS_${VERSION_MUMPS}.tar.gz"
    )
fi

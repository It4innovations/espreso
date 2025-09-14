#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

PREFIX="${1}"
COMPILER_C="${2}"

ELFUTILS_ROOT="${DEPENDENCIES_DIR}/elfutils"
if [ ! -d "${ELFUTILS_ROOT}" ]
then
    sh env/dependencies/clone.elfutils.sh
fi

INSTALL_DIR="${ELFUTILS_ROOT}/install_${PREFIX}"
if [ ! -d "${INSTALL_DIR}" ]
then
    (
        cd "${ELFUTILS_ROOT}"
        mkdir "build_${PREFIX}"
        cd "build_${PREFIX}"
        CC=${COMPILER_C} ../configure --prefix="${INSTALL_DIR}"
        make -j$(nproc)
        make install
    )
fi

export PATH="${INSTALL_DIR}/bin:${PATH}"
export CPATH="${INSTALL_DIR}/include:${CPATH}"
export LIBRARY_PATH="${INSTALL_DIR}/lib:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${INSTALL_DIR}/lib:${LD_LIBRARY_PATH}"

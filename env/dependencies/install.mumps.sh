#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

VERSION_MUMPS=5.8.1

PREFIX="${1}"
COMPILER_C="${2}"
COMPILER_F="${3}"
MY_MUMPS_LAPACK_LINK="${4}"
MY_MUMPS_SCALAP_LINK="${5}"
MY_MUMPS_BLAS_LINK="${6}"

MUMPS_DIR="MUMPS_${VERSION_MUMPS}"
MUMPS_ROOT="${DEPENDENCIES_DIR}/${MUMPS_DIR}"
if [ ! -d "${MUMPS_ROOT}" ]
then
  sh env/dependencies/clone.mumps.sh "${VERSION_MUMPS}"
fi

INSTALL_DIR="${MUMPS_ROOT}/install_${PREFIX}_${COMPILER_C}"
if [ ! -d "${INSTALL_DIR}" ]
then
    (
        ORIG_DIR="${PWD}"
        cd "${MUMPS_ROOT}"
        # originally from Makefile.inc.generic.SEQ
        cp "${ORIG_DIR}/env/dependencies/mumps.Makefile.inc" Makefile.inc
        # makefile uses CC and FC env vars, and also the MY_MUMPS_*_LINK
        export CC="${COMPILER_C}"
        export FC="${COMPILER_F}"
        export MY_MUMPS_LAPACK_LINK
        export MY_MUMPS_SCALAP_LINK
        export MY_MUMPS_BLAS_LINK
        make clean
        make allshared -j$(nproc)
        mkdir -p "${INSTALL_DIR}"
        cd "${INSTALL_DIR}"
        cp -r ../include .
        cp -r ../lib .
    )
fi

prepend_to_CPATH="${INSTALL_DIR}/include"
prepend_to_LIBRARY_PATH="${INSTALL_DIR}/lib"
prepend_to_LD_LIBRARY_PATH="${INSTALL_DIR}/lib"
export CPATH="${prepend_to_CPATH}:${CPATH}"
export LIBRARY_PATH="${prepend_to_LIBRARY_PATH}:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${prepend_to_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}"
echo "          CPATH+=${prepend_to_CPATH}"
echo "   LIBRARY_PATH+=${prepend_to_LIBRARY_PATH}"
echo "LD_LIBRARY_PATH+=${prepend_to_LD_LIBRARY_PATH}"

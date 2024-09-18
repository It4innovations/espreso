#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"

VERSION_SUITESPARSE=v7.6.0

COMPILER_C="${1}"
COMPILER_F="${2}"
MORE_CMAKE_OPTIONS="${3}"

SUITESPARSE_DIR="SuiteSparse_${VERSION_SUITESPARSE}"
SUITESPARSE_ROOT="${DEPENDENCIES_DIR}/${SUITESPARSE_DIR}"
if [ ! -d "${SUITESPARSE_ROOT}" ]
then
  sh ${DEPENDENCIES_DIR}/clone.suitesparse.sh
fi

INSTALL_DIR="${SUITESPARSE_ROOT}/install_${COMPILER_C}"
if [ ! -d "${INSTALL_DIR}" ]
then
    (
        cd "${SUITESPARSE_ROOT}"
        mkdir -p "build_${COMPILER_C}"
        cd "build_${COMPILER_C}"
        cmake -DCMAKE_C_COMPILER="${COMPILER_C}" -DCMAKE_Fortran_COMPILER="${COMPILER_F}" -DENABLE_CUDA=false -DSUITESPARSE_ENABLE_PROJECTS="cholmod;umfpack" ${MORE_CMAKE_OPTIONS} -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} ..
        cmake --build . -j $(nproc)
        cmake --install .
    )
fi

prepend_to_CPATH="${INSTALL_DIR}/include"
prepend_to_LIBRARY_PATH="${INSTALL_DIR}/lib:${INSTALL_DIR}/lib64"
prepend_to_LD_LIBRARY_PATH="${INSTALL_DIR}/lib:${INSTALL_DIR}/lib64"
export CPATH="${prepend_to_CPATH}:${CPATH}"
export LIBRARY_PATH="${prepend_to_LIBRARY_PATH}:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${prepend_to_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}"
echo "          CPATH+=${prepend_to_CPATH}"
echo "   LIBRARY_PATH+=${prepend_to_LIBRARY_PATH}"
echo "LD_LIBRARY_PATH+=${prepend_to_LD_LIBRARY_PATH}"

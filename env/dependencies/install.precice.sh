#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

PREFIX="${1}"
COMPILER_C="${2}"

PRECICE_ROOT="${DEPENDENCIES_DIR}/precice-3.1.2"
if [ ! -d "${PRECICE_ROOT}" ]
then
    sh env/dependencies/clone.precice.sh
fi

INSTALL_DIR="${PRECICE_ROOT}/install_${PREFIX}_${COMPILER_C}"
if [ ! -d "${INSTALL_DIR}" ]
then
    (
        cd "${PRECICE_ROOT}"
        rm -rf build
        mkdir build
        cd build
        cmake .. -DPRECICE_FEATURE_PETSC_MAPPING=OFF -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" -DCMAKE_BUILD_TYPE=Release
        make -j$(nproc)
        make install
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

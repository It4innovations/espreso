#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

PREFIX="${1}"
COMPILER_C="${2}"

GKLIB_ROOT="${DEPENDENCIES_DIR}/gklib"
if [ ! -d "${GKLIB_ROOT}" ]
then
    sh env/dependencies/clone.gklib.sh
fi

INSTALL_DIR="${GKLIB_ROOT}/install_${PREFIX}_${COMPILER_C}"
if [ ! -d "${INSTALL_DIR}" ]
then
    (
        cd "${GKLIB_ROOT}"
        # fix the bug in the compilation script
        sed -i 's/BUILD_SHARED_LIBS/SHARED/'g CMakeLists.txt
        make clean
        make config shared=1 cc=${COMPILER_C} prefix="${INSTALL_DIR}"
        make -j$(nproc)
        make install
        cd "${INSTALL_DIR}/lib64"
        ln -s "libGKlib.so.0" "libGKlib.so"
    )
fi

export CPATH="${INSTALL_DIR}/include:${CPATH}"
export LIBRARY_PATH="${INSTALL_DIR}/lib:${INSTALL_DIR}/lib64:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${INSTALL_DIR}/lib:${INSTALL_DIR}/lib64:${LD_LIBRARY_PATH}"

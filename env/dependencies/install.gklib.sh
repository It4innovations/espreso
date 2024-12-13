#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

COMPILER_C="${1}"

GKLIB_ROOT="${DEPENDENCIES_DIR}/gklib"
if [ ! -d "${GKLIB_ROOT}" ]
then
    sh env/dependencies/clone.gklib.sh
fi

INSTALL_DIR="${GKLIB_ROOT}/install_${COMPILER_C}"
if [ ! -d "${INSTALL_DIR}" ]
then
    (
        cd "${GKLIB_ROOT}"
        sed -i 's/option(BUILD_SHARED_LIBS \"Build shared libraries (.dll\/.so) instead of static ones (.lib\/.a)\" OFF)/option(BUILD_SHARED_LIBS \"Turn on due to intel\" ON)/'g CMakeLists.txt
        make clean
        make config cc=${COMPILER_C} prefix="${INSTALL_DIR}"
        make -j$(nproc)
        make install
    )
fi

export CPATH="${INSTALL_DIR}/include:${CPATH}"
export LIBRARY_PATH="${INSTALL_DIR}/lib:${INSTALL_DIR}/lib64:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${INSTALL_DIR}/lib:${INSTALL_DIR}/lib64:${LD_LIBRARY_PATH}"

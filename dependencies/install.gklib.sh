#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"

GKLIB_ROOT="${DEPENDENCIES_DIR}/gklib"
if [ ! -d "${GKLIB_ROOT}" ]
then
    sh ${DEPENDENCIES_DIR}/clone.gklib.sh
fi

if [ ! -d "${GKLIB_ROOT}/$1" ]
then
    (
        cd "${GKLIB_ROOT}"
        sed -i 's/option(BUILD_SHARED_LIBS \"Build shared libraries (.dll\/.so) instead of static ones (.lib\/.a)\" OFF)/option(BUILD_SHARED_LIBS \"Turn on due to intel\" ON)/'g CMakeLists.txt
        make config cc=$1 prefix="${PWD}/$1"
        make -j$(nproc)
        make install
    )
fi

export CPATH="${GKLIB_ROOT}/$1/include:${CPATH}"
export LIBRARY_PATH="${GKLIB_ROOT}/$1/lib:${GKLIB_ROOT}/$1/lib64:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${GKLIB_ROOT}/$1/lib:${GKLIB_ROOT}/$1/lib64:${LD_LIBRARY_PATH}"

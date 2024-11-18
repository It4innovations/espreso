#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"

METIS_ROOT="${DEPENDENCIES_DIR}/metis"
if [ ! -d "${METIS_ROOT}" ]
then
    sh ${DEPENDENCIES_DIR}/clone.metis.sh
fi

if [ ! -d "${METIS_ROOT}/$1_32" ]
then
    (
        cd "${METIS_ROOT}"
        sed -i 's/add_subdirectory(\"programs\")/#add_subdirectory(\"programs\")/'g CMakeLists.txt
        make config shared=1 cc=$1 prefix="${PWD}/$1_32"
        make -j$(nproc)
        make install
    )
fi

export CPATH="${METIS_ROOT}/$1_32/include:${CPATH}"
export LIBRARY_PATH="${METIS_ROOT}/$1_32/lib:${METIS_ROOT}/$1_32/lib64:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${METIS_ROOT}/$1_32/lib:${METIS_ROOT}/$1_32/lib64:${LD_LIBRARY_PATH}"
echo "          CPATH+=${METIS_ROOT}/$1_32/include"
echo "   LIBRARY_PATH+=${METIS_ROOT}/$1_32/lib"
echo "LD_LIBRARY_PATH+=${METIS_ROOT}/$1_32/lib"

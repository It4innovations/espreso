#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"

PRECICE_ROOT="${DEPENDENCIES_DIR}/precice-3.1.2"
if [ ! -d "${PRECICE_ROOT}" ]
then
    sh ${DEPENDENCIES_DIR}/clone.precice.sh
fi

if [ ! -d "${PRECICE_ROOT}/$1" ]
then
    (
        cd "${PRECICE_ROOT}"
        rm -rf build
        mkdir build
        cd build
        cmake .. -DPRECICE_FEATURE_PETSC_MAPPING=OFF -DCMAKE_INSTALL_PREFIX="${PRECICE_ROOT}/$1" -DCMAKE_BUILD_TYPE=Release
        make -j$(nproc)
        make install
    )
fi

export CPATH="${PRECICE_ROOT}/$1/include:${CPATH}"
export LIBRARY_PATH="${PRECICE_ROOT}/$1/lib:${PRECICE_ROOT}/$1/lib64:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${PRECICE_ROOT}/$1/lib:${PRECICE_ROOT}/$1/lib64:${LD_LIBRARY_PATH}"
echo "          CPATH+=${PRECICE_ROOT}/$1/include"
echo "   LIBRARY_PATH+=${PRECICE_ROOT}/$1/lib"
echo "LD_LIBRARY_PATH+=${PRECICE_ROOT}/$1/lib"

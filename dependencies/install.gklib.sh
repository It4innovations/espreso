#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"

GKLIB_ROOT="${DEPENDENCIES_DIR}/gklib"
if [ ! -d "${GKLIB_ROOT}" ]
then
    sh ${DEPENDENCIES_DIR}/clone.gklib.sh
fi

if [ ! -d "${GKLIB_ROOT}/include" ]
then
    (
        cd "${GKLIB_ROOT}"
        make config shared=1 cc=$1 prefix="${PWD}"
        make -j$(nproc)
        make install
    )
fi

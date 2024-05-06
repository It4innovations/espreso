#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"

PARMETIS_ROOT="${DEPENDENCIES_DIR}/parmetis-4.0.3"
if [ ! -d "${PARMETIS_ROOT}" ]
then
    sh ${DEPENDENCIES_DIR}/clone.parmetis.sh
fi

if [ ! -d "${PARMETIS_ROOT}/$1_32" ]
then
    (
        cd "${PARMETIS_ROOT}"
        make config shared=1 cc=$1 prefix="${PWD}/$1_32"
        make -j$(nproc)
        make install
    )
fi

export CPATH="${PARMETIS_ROOT}/$1_32/include:${CPATH}"
export LIBRARY_PATH="${PARMETIS_ROOT}/$1_32/lib:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${PARMETIS_ROOT}/$1_32/lib:${LD_LIBRARY_PATH}"

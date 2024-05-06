#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"

METIS_ROOT="${DEPENDENCIES_DIR}/metis-5.1.0"
if [ ! -d "${METIS_ROOT}" ]
then
    sh ${DEPENDENCIES_DIR}/clone.metis.sh
fi

if [ ! -d "${METIS_ROOT}/$1_64" ]
then
    (
        cd "${METIS_ROOT}"
        sed -i 's/#define IDXTYPEWIDTH 32/#define IDXTYPEWIDTH 64/g' ${PWD}/include/metis.h

        make config shared=1 cc=$1 gklib_path="${PWD}/GKlib" prefix="${PWD}/$1_64"
        make -j$(nproc)
        make install
    )
fi

export CPATH="${METIS_ROOT}/$1_64/include:${CPATH}"
export LIBRARY_PATH="${METIS_ROOT}/$1_64/lib:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${METIS_ROOT}/$1_64/lib:${LD_LIBRARY_PATH}"
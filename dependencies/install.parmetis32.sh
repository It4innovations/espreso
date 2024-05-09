#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"

METIS_ROOT="${DEPENDENCIES_DIR}/metis-5.1.0"
PARMETIS_ROOT="${DEPENDENCIES_DIR}/parmetis-4.0.3"
if [ ! -d "${PARMETIS_ROOT}" ]
then
    sh ${DEPENDENCIES_DIR}/clone.parmetis.sh
fi

if [ ! -d "${PARMETIS_ROOT}/$1_32" ]
then
    (
        cd "${PARMETIS_ROOT}"
        # skip building of METIS
        sed -i 's/add_subdirectory(${METIS_PATH}/#add_subdirectory(${METIS_PATH}/'g CMakeLists.txt
        sed -i 's/add_subdirectory(programs)/#add_subdirectory(programs)/'g CMakeLists.txt
        make config shared=1 cc=$1 prefix="${PWD}/$1_32" metis_path=${METIS_ROOT} gklib_path=${METIS_ROOT}/GKlib
        make -j$(nproc)
        make install
    )
fi

export CPATH="${PARMETIS_ROOT}/$1_32/include:${CPATH}"
export LIBRARY_PATH="${PARMETIS_ROOT}/$1_32/lib:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${PARMETIS_ROOT}/$1_32/lib:${LD_LIBRARY_PATH}"

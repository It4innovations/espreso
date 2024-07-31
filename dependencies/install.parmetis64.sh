#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"

PARMETIS_ROOT="${DEPENDENCIES_DIR}/parmetis-4.0.3"
if [ ! -d "${PARMETIS_ROOT}" ]
then
    sh ${DEPENDENCIES_DIR}/clone.parmetis.sh
fi

if [ ! -d "${PARMETIS_ROOT}/$1_64" ]
then
    (
        cd "${PARMETIS_ROOT}"
        sed -i 's/add_subdirectory(${METIS_PATH}/#add_subdirectory(${METIS_PATH}/'g CMakeLists.txt
        sed -i 's/add_subdirectory(programs)/#add_subdirectory(programs)/'g CMakeLists.txt
        sed -i 's/#define IDXTYPEWIDTH 32/#define IDXTYPEWIDTH 64/g' ${PARMETIS_ROOT}/metis/include/metis.h
        make config shared=1 cc=$1 gklib_path="${DEPENDENCIES_DIR}/gklib" prefix="${PWD}/$1_64"
        make -j$(nproc)
        make install
    )
fi

export CPATH="${PARMETIS_ROOT}/$1_64/include:${CPATH}"
export LIBRARY_PATH="${PARMETIS_ROOT}/$1_64/lib:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${PARMETIS_ROOT}/$1_64/lib:${LD_LIBRARY_PATH}"
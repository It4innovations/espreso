#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"

PARMETIS_ROOT="${DEPENDENCIES_DIR}/parmetis"
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
        make config shared=1 cc=$1 prefix="${PWD}/$1_64"
        make -j$(nproc)
        make install
    )
fi

export CPATH="${PARMETIS_ROOT}/$1_64/include:${CPATH}"
export LIBRARY_PATH="${PARMETIS_ROOT}/$1_64/lib:${PARMETIS_ROOT}/$1_32/lib64:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${PARMETIS_ROOT}/$1_64/lib:${PARMETIS_ROOT}/$1_32/lib64:${LD_LIBRARY_PATH}"
echo "          CPATH+=${PARMETIS_ROOT}/$1_64/include"
echo "   LIBRARY_PATH+=${PARMETIS_ROOT}/$1_64/lib"
echo "LD_LIBRARY_PATH+=${PARMETIS_ROOT}/$1_64/lib"

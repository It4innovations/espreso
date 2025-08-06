#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

PREFIX="${1}"
COMPILER_C="${2}"

PARMETIS_ROOT="${DEPENDENCIES_DIR}/parmetis"
if [ ! -d "${PARMETIS_ROOT}" ]
then
    sh env/dependencies/clone.parmetis.sh
fi

INSTALL_DIR="${PARMETIS_ROOT}/install_${PREFIX}_${COMPILER_C}_64"
if [ ! -d "${INSTALL_DIR}" ]
then
    (
        cd "${PARMETIS_ROOT}"
        sed -i 's/add_subdirectory(${METIS_PATH}/#add_subdirectory(${METIS_PATH}/'g CMakeLists.txt
        sed -i 's/add_subdirectory(programs)/#add_subdirectory(programs)/'g CMakeLists.txt
        make clean
        make config shared=1 cc=${COMPILER_C} prefix="${INSTALL_DIR}"
        make -j$(nproc)
        make install
    )
fi

prepend_to_CPATH="${INSTALL_DIR}/include"
prepend_to_LIBRARY_PATH="${INSTALL_DIR}/lib:${INSTALL_DIR}/lib64"
prepend_to_LD_LIBRARY_PATH="${INSTALL_DIR}/lib:${INSTALL_DIR}/lib64"
export CPATH="${prepend_to_CPATH}:${CPATH}"
export LIBRARY_PATH="${prepend_to_LIBRARY_PATH}:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${prepend_to_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}"
echo "          CPATH+=${prepend_to_CPATH}"
echo "   LIBRARY_PATH+=${prepend_to_LIBRARY_PATH}"
echo "LD_LIBRARY_PATH+=${prepend_to_LD_LIBRARY_PATH}"


#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

VERSION_PASTIX=v6.4.0

PREFIX="${1}"
COMPILER_CXX="${2}"
COMPILER_C="${3}"
COMPILER_FORTRAN="${4}"
BLA_VENDOR="${5}"

PASTIX_DIR="pastix_${VERSION_PASTIX}"
PASTIX_ROOT="${DEPENDENCIES_DIR}/${PASTIX_DIR}"
if [ ! -d "${PASTIX_ROOT}" ]
then
    sh env/dependencies/clone.pastix.sh
fi



INSTALL_DIR="${PASTIX_ROOT}/install_${PREFIX}_${COMPILER_C}"
if [ ! -d "${INSTALL_DIR}" ]
then
    (
        cd "${PASTIX_ROOT}"
        mkdir -p "build_${PREFIX}_${COMPILER_C}"
        cd "build_${PREFIX}_${COMPILER_C}"
        cmake -DBLA_VENDOR="${BLA_VENDOR}" -DMETIS_LIBRARIES="-lmetis -lGKlib -lm" -DPASTIX_WITH_MPI=OFF -DPASTIX_WITH_FORTRAN=OFF -DPASTIX_INT64=OFF -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" -DCMAKE_C_COMPILER="${COMPILER_C}" -DCMAKE_CXX_COMPILER="${COMPILER_CXX}" -DCMAKE_Fortran_COMPILER="${COMPILER_FORTRAN}" -DPASTIX_WITH_CUDA=OFF -DPASTIX_WITH_PARSEC=OFF -DPASTIX_WITH_STARPU=OFF -DPASTIX_ORDERING_METIS=ON -DPASTIX_ORDERING_SCOTCH=OFF -DBUILD_SHARED_LIBS=ON ..
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


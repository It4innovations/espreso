#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

VERSION_SUPERLU=v9.1.0

PREFIX="${1}"
COMPILER_CXX="${2}"
COMPILER_C="${3}"
MY_SUPERLU_BLAS_LINK="${4}"

SUPERLU_DIR="SuperLU_dist_${VERSION_SUPERLU}"
SUPERLU_ROOT="${DEPENDENCIES_DIR}/${SUPERLU_DIR}"
if [ ! -d "${SUPERLU_ROOT}" ]
then
  sh env/dependencies/clone.superlu_dist.sh
fi

INSTALL_DIR="${SUPERLU_ROOT}/install_${PREFIX}_${COMPILER_C}"
if [ ! -d "${INSTALL_DIR}" ]
then
    (
        cd "${SUPERLU_ROOT}"
        mkdir -p "build_${PREFIX}_${COMPILER_C}"
        cd "build_${PREFIX}_${COMPILER_C}"
        cmake -DTPL_BLAS_LIBRARIES="${MY_SUPERLU_BLAS_LINK}" -DTPL_ENABLE_LAPACKLIB=ON -DTPL_ENABLE_CUDALIB=TRUE -DCMAKE_CUDA_ARCHITECTURES=native -DTPL_ENABLE_PARMETISLIB=ON -DTPL_PARMETIS_INCLUDE_DIRS="." -DTPL_PARMETIS_LIBRARIES="-lparmetis -lmetis" -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" -DCMAKE_C_COMPILER="${COMPILER_C}" -DCMAKE_CXX_COMPILER="${COMPILER_CXX}" -DBUILD_SHARED_LIBS=ON ..
        make -j$(nproc)
        make install
    )
fi

prepend_to_CPATH="${INSTALL_DIR}/include"
prepend_to_LIBRARY_PATH="${INSTALL_DIR}/lib64"
prepend_to_LD_LIBRARY_PATH="${INSTALL_DIR}/lib64"
export CPATH="${prepend_to_CPATH}:${CPATH}"
export LIBRARY_PATH="${prepend_to_LIBRARY_PATH}:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${prepend_to_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}"
echo "          CPATH+=${prepend_to_CPATH}"
echo "   LIBRARY_PATH+=${prepend_to_LIBRARY_PATH}"
echo "LD_LIBRARY_PATH+=${prepend_to_LD_LIBRARY_PATH}"

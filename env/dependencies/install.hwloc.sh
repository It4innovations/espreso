#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

PREFIX="${1}"
COMPILER_C="${2}"

HWLOC_DIR="hwloc"
HWLOC_ROOT="${DEPENDENCIES_DIR}/${HWLOC_DIR}"
if [ ! -d "${HWLOC_ROOT}" ]
then
    sh env/dependencies/clone.hwloc.sh
fi



INSTALL_DIR="${HWLOC_ROOT}/install_${PREFIX}_${COMPILER_C}"
if [ ! -d "${INSTALL_DIR}" ]
then
    (
        cd "${HWLOC_ROOT}"
        mkdir "build_${PREFIX}_${COMPILER_C}"
        cd "build_${PREFIX}_${COMPILER_C}"
        ../autogen.sh
        CC="${COMPILER_C}" ../configure --prefix="${INSTALL_DIR}"
        make -j$(nproc)
        make install
    )
fi



prepend_to_PATH="${INSTALL_DIR}/bin"
prepend_to_CPATH="${INSTALL_DIR}/include"
prepend_to_LIBRARY_PATH="${INSTALL_DIR}/lib"
prepend_to_LD_LIBRARY_PATH="${INSTALL_DIR}/lib"
prepend_to_PKG_CONFIG_PATH="${INSTALL_DIR}/lib/pkgconfig"
export PATH="${prepend_to_PATH}:${PATH}"
export CPATH="${prepend_to_CPATH}:${CPATH}"
export LIBRARY_PATH="${prepend_to_LIBRARY_PATH}:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${prepend_to_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}"
export PKG_CONFIG_PATH="${prepend_to_PKG_CONFIG_PATH}:${PKG_CONFIG_PATH}"
echo "           PATH+=${prepend_to_PATH}"
echo "          CPATH+=${prepend_to_CPATH}"
echo "   LIBRARY_PATH+=${prepend_to_LIBRARY_PATH}"
echo "LD_LIBRARY_PATH+=${prepend_to_LD_LIBRARY_PATH}"
echo "PKG_CONFIG_PATH+=${prepend_to_PKG_CONFIG_PATH}"

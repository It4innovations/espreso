#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

AOCL_ROOT="${DEPENDENCIES_DIR}/aocl-linux-aocc-4.1.0"
if [ ! -d "${AOCL_ROOT}" ]
then
    sh env/dependencies/clone.aocl.sh
fi

INSTALL_DIR="${AOCL_ROOT}/install"
if [ ! -d "${INSTALL_DIR}" ]
then
    (
        cd "${AOCL_ROOT}"
        ./install.sh -t "${INSTALL_DIR}" -i lp64
    )
fi
source "${INSTALL_DIR}/4.1.0/aocc/amd-libs.cfg"

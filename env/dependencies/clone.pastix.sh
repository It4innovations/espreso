#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

VERSION_PASTIX=v6.4.0

PASTIX_DIR="pastix_${VERSION_PASTIX}"
PASTIX_ROOT="${DEPENDENCIES_DIR}/${PASTIX_DIR}"
if [ ! -d "${PASTIX_ROOT}" ]
then
    (
        echo "pastix ${VERSION_PASTIX} not found, cloning..."
        cd "${DEPENDENCIES_DIR}"
        git clone -b "${VERSION_PASTIX}" https://gitlab.inria.fr/solverstack/pastix.git "${PASTIX_DIR}"
        cd "${PASTIX_ROOT}"
        git submodule init
        git submodule update
        mv .git .git_ # so that git submodule update is not executed again in install cmake. good for offline compilation.
    )
fi

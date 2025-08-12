#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

VERSION_PASTIX=release-6.4.0

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
    )
fi

# pastix requires this specific branch of parsec
PARSEC_DIR="pastix_parsec"
PARSEC_ROOT="${PASTIX_ROOT}/${PARSEC_DIR}"
if [ ! -d "${PARSEC_ROOT}" ]
then
    (
        cd "${PASTIX_ROOT}"
        git clone -b pastix-6.0.2 https://bitbucket.org/mfaverge/parsec "${PARSEC_DIR}"
        cd "${PARSEC_ROOT}"
    )
fi

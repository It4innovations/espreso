#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

AOCL_ROOT="${DEPENDENCIES_DIR}/aocl-linux-aocc-4.1.0"
if [ ! -d "${AOCL_ROOT}" ]
then
    (
        echo "AOCL 4.1.0 not found, downloading..."
        cd "${DEPENDENCIES_DIR}"
        ../env/dependencies/download.sh https://download.amd.com/developer/eula/aocl/aocl-4-1/aocl-linux-aocc-4.1.0.tar.gz
        tar -xf aocl-linux-aocc-4.1.0.tar.gz
        rm aocl-linux-aocc-4.1.0.tar.gz
    )
fi

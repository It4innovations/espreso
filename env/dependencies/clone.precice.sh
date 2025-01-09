#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

PRECICE_ROOT="${DEPENDENCIES_DIR}/precice-3.1.2"
if [ ! -d "${PRECICE_ROOT}" ]
then
    (
        echo "PreCICE not installed, cloning..."
        cd "${DEPENDENCIES_DIR}"
        ../env/dependencies/download.sh https://github.com/precice/precice/archive/v3.1.2.tar.gz
        tar -xzvf v3.1.2.tar.gz
        rm v3.1.2.tar.gz
    )
fi
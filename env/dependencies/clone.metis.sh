#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

METIS_ROOT="${DEPENDENCIES_DIR}/metis"
if [ ! -d "${METIS_ROOT}" ]
then
    (
        echo "METIS not installed, cloning..."
        cd "${DEPENDENCIES_DIR}"
        git clone https://github.com/KarypisLab/METIS.git metis
        cd metis
        git checkout e0f1b88b8efcb24ffa0ec55eabb78fbe61e58ae7
    )
fi

#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

PARMETIS_ROOT="${DEPENDENCIES_DIR}/parmetis"
if [ ! -d "${PARMETIS_ROOT}" ]
then
    (
        echo "ParMETIS not installed, cloning..."
        cd "${DEPENDENCIES_DIR}"
        git clone https://github.com/KarypisLab/ParMETIS.git parmetis
    )
fi

#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"

PARMETIS_ROOT="${DEPENDENCIES_DIR}/parmetis-4.0.3"
if [ ! -d "${PARMETIS_ROOT}" ]
then
    (
        echo "ParMETIS not installed, cloning..."
        cd "${DEPENDENCIES_DIR}"
        git clone https://github.com/KarypisLab/ParMETIS.git parmetis
    )
fi
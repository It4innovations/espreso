#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"

GKLIB_ROOT="${DEPENDENCIES_DIR}/gklib"
if [ ! -d "${GKLIB_ROOT}" ]
then
    (
        echo "GKlib not installed, cloning..."
        cd "${DEPENDENCIES_DIR}"
        git clone https://github.com/KarypisLab/GKlib.git gklib
    )
fi
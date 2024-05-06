#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"

METIS_ROOT="${DEPENDENCIES_DIR}/metis-5.1.0"
if [ ! -d "${METIS_ROOT}" ]
then
    (
        echo "METIS not installed, cloning..."
        cd "${DEPENDENCIES_DIR}"
        curl -s -O http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
        tar -xf metis-5.1.0.tar.gz
        rm metis-5.1.0.tar.gz
    )
fi

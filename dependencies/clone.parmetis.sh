#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"

PARMETIS_ROOT="${DEPENDENCIES_DIR}/parmetis-4.0.3"
if [ ! -d "${PARMETIS_ROOT}" ]
then
    (
        echo "ParMETIS not installed, cloning..."
        cd "${DEPENDENCIES_DIR}"
        curl -s -O http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
        tar -xf parmetis-4.0.3.tar.gz
        rm parmetis-4.0.3.tar.gz
    )
fi
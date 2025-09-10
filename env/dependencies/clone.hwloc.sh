#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

HWLOC_ROOT="${DEPENDENCIES_DIR}/hwloc"
if [ ! -d "${HWLOC_ROOT}" ]
then
    (
        echo "hwloc not installed, cloning..."
        cd "${DEPENDENCIES_DIR}"
        git clone -b hwloc-2.12.2 https://github.com/open-mpi/hwloc.git hwloc
    )
fi

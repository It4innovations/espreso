#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"

version_name="nvpl-linux-sbsa-24.7"

NVPL_DIR="${version_name}"
NVPL_ROOT="${DEPENDENCIES_DIR}/${NVPL_DIR}"
if [ ! -d "${NVPL_ROOT}" ]
then
    (
        echo "NVPL not found, downloading..."
        cd "${DEPENDENCIES_DIR}"
        wget https://developer.download.nvidia.com/compute/nvpl/24.7/local_installers/${version_name}.tar.gz
        tar -xf "${version_name}.tar.gz"
        rm "${version_name}.tar.gz"
        cd "${version_name}/include"
        ln -s nvpl_blas_cblas.h cblas.h
    )
fi


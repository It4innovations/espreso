#!/bin/bash

if [ $# -lt 1 ]
then
    echo "architecture not specified"
    return 1
fi
arch="${1}"

DEPENDENCIES_DIR="${PWD}/dependencies"

VERSION_CMAKE="3.30.5"
version_name="cmake-${VERSION_CMAKE}-linux-${arch}"
tarball_url="https://github.com/Kitware/CMake/releases/download/v${VERSION_CMAKE}/${version_name}.tar.gz"

CMAKE_DIR="${version_name}"
CMAKE_ROOT="${DEPENDENCIES_DIR}/${CMAKE_DIR}"
if [ ! -d "${CMAKE_ROOT}" ]
then
    (
        echo "CMake ${VERSION_CMAKE} not found, downloading..."
        cd "${DEPENDENCIES_DIR}"
        wget "${tarball_url}"
        tar -xf "${version_name}.tar.gz"
        rm "${version_name}.tar.gz"
    )
fi


#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"

VERSION_CMAKE="3.30.5"
version_name="cmake-${VERSION_CMAKE}-linux-aarch64"

CMAKE_DIR="${version_name}"
CMAKE_ROOT="${DEPENDENCIES_DIR}/${CMAKE_DIR}"
if [ ! -d "${CMAKE_ROOT}" ]
then
  sh ${DEPENDENCIES_DIR}/clone.cmake.sh
fi

export PATH="${CMAKE_ROOT}/bin:${PATH}"

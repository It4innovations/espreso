#!/bin/bash

if [ $# -lt 1 ]
then
  echo "architecture not specified"
  return 1
fi
arch="${1}"

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

VERSION_CMAKE="3.30.5"
version_name="cmake-${VERSION_CMAKE}-linux-${arch}"

CMAKE_DIR="${version_name}"
CMAKE_ROOT="${DEPENDENCIES_DIR}/${CMAKE_DIR}"
if [ ! -d "${CMAKE_ROOT}" ]
then
  sh ${DEPENDENCIES_DIR}/clone.cmake.sh "${arch}"
fi

export PATH="${CMAKE_ROOT}/bin:${PATH}"

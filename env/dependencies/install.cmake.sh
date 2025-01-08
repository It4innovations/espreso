#!/bin/bash

if [ $# -lt 2 ]
then
  echo "not enough arguments"
  return 1
fi
prefix="${1}"
arch="${2}"

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

VERSION_CMAKE="3.30.5"
version_name="cmake-${VERSION_CMAKE}-linux-${arch}"

CMAKE_DIR="${version_name}"
CMAKE_ROOT="${DEPENDENCIES_DIR}/${CMAKE_DIR}"
if [ ! -d "${CMAKE_ROOT}" ]
then
  sh env/dependencies/clone.cmake.sh "${arch}"
fi

export PATH="${CMAKE_ROOT}/bin:${PATH}"

#!/bin/bash

PREFIX="${1}"

DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

NOSE_DIR="${DEPENDENCIES_DIR}/nose2/${PREFIX}"
if [ ! -d "${NOSE_DIR}" ]
then
  pip install nose2 --target=${NOSE_DIR}
fi

export PATH="${NOSE_DIR}/bin:${PATH}"

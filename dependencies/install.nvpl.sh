#!/bin/bash

DEPENDENCIES_DIR="${PWD}/dependencies"

version_name="nvpl-linux-sbsa-24.7"

NVPL_DIR="${version_name}"
NVPL_ROOT="${DEPENDENCIES_DIR}/${NVPL_DIR}"
if [ ! -d "${NVPL_ROOT}" ]
then
  sh ${DEPENDENCIES_DIR}/clone.nvpl.sh
fi

export CPATH="${NVPL_ROOT}/include:${CPATH}"
export LIBRARY_PATH="${NVPL_ROOT}/lib:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${NVPL_ROOT}/lib:${LD_LIBRARY_PATH}"

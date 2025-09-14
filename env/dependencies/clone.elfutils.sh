#!/bin/bash

ESP_ROOT_DIR="${PWD}"
DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

ELFUTILS_ROOT="${DEPENDENCIES_DIR}/elfutils"
if [ ! -d "${ELFUTILS_ROOT}" ]
then
    (
        echo "elfutils not installed, cloning..."
        cd "${DEPENDENCIES_DIR}"
        wget https://sourceware.org/elfutils/ftp/0.193/elfutils-0.193.tar.bz2
        tar -xf elfutils-0.193.tar.bz2
        rm elfutils-0.193.tar.bz2
        mv elfutils-0.193 elfutils
    )
fi
#!/bin/bash

ORIG_DIR="${PWD}"
DEPENDENCIES_DIR="${PWD}/dependencies"
mkdir -p "${DEPENDENCIES_DIR}"

GKLIB_ROOT="${DEPENDENCIES_DIR}/gklib"
if [ ! -d "${GKLIB_ROOT}" ]
then
    (
        echo "GKlib not installed, cloning..."
        cd "${DEPENDENCIES_DIR}"
        git clone https://github.com/KarypisLab/GKlib.git gklib
        cd gklib
        checkout 804361cce8dbf2459863d329498dbe0f2d34d68e
        ln -s "${ORIG_DIR}/env/dependencies/patch_gklib.txt"
        perl -i -p0e 's/\/\* The array for the state vector \*\/\nstatic uint64_t mt\[NN\]; \n\/\* mti==NN\+1 means mt\[NN\] is not initialized \*\/\nstatic int mti=NN\+1; \n/`cat patch_gklib.txt`/se' src/random.c
        sed -i '1i#define USE_GKRAND' src/random.c
        rm patch_gklib.txt
    )
fi
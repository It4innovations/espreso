#!/bin/bash

if [[ "$(which icpx)" != *"2024.2"* ]]
then
    echo
    echo "Warning, default Intel toolkit changed, things might work differently"
    echo
fi



. dependencies/install.suitesparse.sh icx
. dependencies/install.gklib.sh icx
. dependencies/install.metis32.sh icx
. dependencies/install.parmetis32.sh mpiicx



# ugly hack, no idea why or how it works, but solves the relink..with error
if [ ! -e build/libirc.so ]
then
    mkdir -p build
    ln -s /lib/x86_64-linux-gnu/libc.so.6 build/libirc.so
fi
export LIBRARY_PATH="${PWD}/build:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${PWD}/build:${LD_LIBRARY_PATH}"
# export LD_PRELOAD="${PWD}/build/libirc.so:${LD_PRELOAD}"



export CXX=mpiicpx
export CXXFLAGS+=" -ferror-limit=1"



export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=12,1

export ESPRESO_FORBID_MKL_PARDISO="1"
export ESPRESO_RANK_TO_GPU_MAP="0"
# export ESPRESO_SYCL_TARGETS="intel_gpu_pvc"
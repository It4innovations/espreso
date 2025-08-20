#!/bin/bash

if [ $# -lt 1 ]
then
    echo "ERROR: You have to specify which cpu spsolver to use (mkl/suitesparse)"
    return 1
fi
cpu_spsolver="${1}"
if [ "${cpu_spsolver}" == "ss" ]; then cpu_spsolver="suitesparse"; fi
if [ "${cpu_spsolver}" == "mklpardiso" ]; then cpu_spsolver="mkl"; fi
if [ "${cpu_spsolver}" != "suitesparse" ] && [ "${cpu_spsolver}" != "mkl" ]
then
    echo "ERROR: wrong cpu spsolver"
    return 2
fi



buildname="build-it4i-cs-spr-${cpu_spsolver}"
export WAFLOCK=".lock-waf_linux_${buildname}"
curr_buildname="$(readlink build)"
if [ "${curr_buildname}" != "${buildname}" ]
then
    rm -rf build
    ln -s "${buildname}" build
fi



if [ ! -d "dependencies/oneAPI-2024.2.0/install" ]
then
    echo "intel toolkit not installed in dependencies, please install it manually"
    return 3
fi
source dependencies/oneAPI-2024.2.0/install/setvars.sh



. dependencies/install.cmake.sh x86_64
if [ "${cpu_spsolver}" == "suitesparse" ]
then
    . dependencies/install.suitesparse.sh icx ifx
fi
. env/dependencies/install.gklib.sh p10 icx
. env/dependencies/install.metis32.sh p10 icx
. env/dependencies/install.parmetis32.sh p10 mpiicx



export CXX=mpiicpx
export CXXFLAGS+=" -ferror-limit=1"

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=12,1

export ESPRESO_USE_WRAPPER_DNBLAS=mkl
export ESPRESO_USE_WRAPPER_DNSOLVER=mkl
export ESPRESO_USE_WRAPPER_LAPACK=mkl
export ESPRESO_USE_WRAPPER_SPBLAS=mkl
export ESPRESO_USE_WRAPPER_SPSOLVER="${cpu_spsolver}"
export ESPRESO_USE_WRAPPER_SCSOLVER="${cpu_spsolver}"

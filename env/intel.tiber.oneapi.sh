#!/bin/bash

if [ $# -lt 1 ]
then
    echo "ERROR: You have to specify which cpu spsolver to use (mklpardiso/suitesparse)"
    return 1
fi
cpu_spsolver="${1}"
if [ "${cpu_spsolver}" == "ss" ]; then cpu_spsolver="suitesparse"; fi
if [ "${cpu_spsolver}" == "mkl" ]; then cpu_spsolver="mklpardiso"; fi
if [ "${cpu_spsolver}" != "suitesparse" ] && [ "${cpu_spsolver}" != "mklpardiso" ]
then
    echo "ERROR: wrong cpu spsolver"
    return 2
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
. dependencies/install.gklib.sh icx
. dependencies/install.metis32.sh icx
. dependencies/install.parmetis32.sh mpiicx



export CXX=mpiicpx
export CXXFLAGS+=" -ferror-limit=1"

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=6,1

if [ "${cpu_spsolver}" == "suitesparse" ]
then
    export ESPRESO_FORBID_MKL_PARDISO="1"
fi
export ESPRESO_RANK_TO_GPU_MAP="0,1,2,3,4,5,6,7"
# in benchmarks I will set the visible devices to only one for each hq worker

export ZE_FLAT_DEVICE_HIERARCHY="FLAT"
# export ESPRESO_SYCL_TARGETS="intel_gpu_pvc"

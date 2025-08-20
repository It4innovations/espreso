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



buildname="build-intel-tiber-oneapi-${cpu_spsolver}"
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



. env/dependencies/install.cmake.sh oneapi x86_64
if [ "${cpu_spsolver}" == "suitesparse" ]
then
    . env/dependencies/install.suitesparse.sh oneapi icx ifx
fi
. env/dependencies/install.gklib.sh oneapi icx
. env/dependencies/install.metis32.sh oneapi icx
. env/dependencies/install.parmetis32.sh oneapi mpiicx



export CXX=mpiicpx
export CXXFLAGS+=" -ferror-limit=1"

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=6,1

# export ESPRESO_RANK_TO_GPU_MAP="$(seq -s, 0 15)"
export ESPRESO_RANK_TO_GPU_MAP="0"
# in benchmarks I will set the visible devices to only one for each hq worker

export ZE_FLAT_DEVICE_HIERARCHY="FLAT" # two stacks on a GPU card appear each as a single GPU device
export ZES_ENABLE_SYSMAN="1" # to make free_memory device info available

#https://www.intel.com/content/www/us/en/docs/oneapi/optimization-guide-gpu/2024-2/ahead-of-time-compilation.html
# export ESPRESO_SYCL_TARGETS="spir64_gen"
# export ESPRESO_SYCL_BACKEND_OPTIONS="-device pvc"
# icpx: warning: argument unused during compilation: '-Xs -device pvc' [-Wunused-command-line-argument]
# i dont understand ...

export ESPRESO_USE_WRAPPER_DNBLAS=mkl
export ESPRESO_USE_WRAPPER_DNSOLVER=mkl
export ESPRESO_USE_WRAPPER_LAPACK=mkl
export ESPRESO_USE_WRAPPER_SPBLAS=mkl
export ESPRESO_USE_WRAPPER_SPSOLVER="${cpu_spsolver}"
export ESPRESO_USE_WRAPPER_SCSOLVER="${cpu_spsolver}"
export ESPRESO_USE_WRAPPER_GPU=oneapi

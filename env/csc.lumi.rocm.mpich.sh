#!/bin/bash

buildname="build-csc-lumi-rocm-mpich"
export WAFLOCK=".lock-waf_linux_${buildname}"
curr_buildname="$(readlink build)"
if [ "${curr_buildname}" != "${buildname}" ]
then
    rm -rf build
    ln -s "${buildname}" build
fi



module load LUMI/24.03
module load partition/G
module load rocm/6.0.3
module load aocc/4.1.0



. env/dependencies/install.cmake.sh rocm x86_64
. env/dependencies/install.aocl.sh
. env/dependencies/install.suitesparse.sh rocm hipcc hipfc "-DLAPACK_LIBRARIES=${PWD}/dependencies/aocl-linux-aocc-4.1.0/install/4.1.0/aocc/lib/libflame.so.4.1.0"
. env/dependencies/install.gklib.sh rocm hipcc
. env/dependencies/install.metis32.sh rocm hipcc
. env/dependencies/install.parmetis32.sh rocm mpicc



export LIBRARIES=mpich
export CPATH="${CRAY_MPICH_DIR}/include:${CPATH}"
export LIBRARY_PATH="${CRAY_MPICH_DIR}/lib:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${CRAY_MPICH_DIR}/lib:${LD_LIBRARY_PATH}"

export CXX=hipcc
export CXXFLAGS+=" -ferror-limit=1"
export BLAS_LIBRARIES=blis
export LAPACK_LIBRARIES=flame

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=7,1

export ESPRESO_RANK_TO_GPU_MAP="0"
# export ESPRESO_RANK_TO_GPU_MAP="4,5,2,3,6,7,0,1" # dont use

export ESPRESO_ROCM_ARCH="gfx90a:sramecc+:xnack-"

export ESPRESO_USE_WRAPPER_DNBLAS=blas
export ESPRESO_USE_WRAPPER_DNSOLVER=lapack
export ESPRESO_USE_WRAPPER_LAPACK=lapack
export ESPRESO_USE_WRAPPER_SPBLAS=suitesparse
export ESPRESO_USE_WRAPPER_SPSOLVER=suitesparse
export ESPRESO_USE_WRAPPER_SCSOLVER=suitesparse
export ESPRESO_USE_WRAPPER_GPU=rocm



# always use slurm options --gpus-per-task=1 --cpus-per-task=7
# salloc --account=project_465000572 --partition=standard-g --ntasks=8 --gpus-per-task=1 --cpus-per-task=7  --time=8:00:00
# srun -n 8 ./build/espreso -c path/to/espreso.ecf

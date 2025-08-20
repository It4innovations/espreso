#!/bin/bash

buildname="build-csc-lumi-rocm-mkl"
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

if [ ! -d "dependencies/oneAPI-2024.2.0/install" ]
then
    echo "intel toolkit not installed in dependencies, please install it manually"
    return 1
fi
source dependencies/oneAPI-2024.2.0/install/setvars.sh



. env/dependencies/install.gklib.sh rocm hipcc
. env/dependencies/install.metis32.sh rocm hipcc
. env/dependencies/install.parmetis32.sh rocm mpicc



export LD_LIBRARY_PATH="/opt/rocm-6.0.3/lib/llvm/bin/../lib:${LD_LIBRARY_PATH}" # workaround

export LIBRARIES=mpich
export CPATH="${CRAY_MPICH_DIR}/include:${CPATH}"
export LIBRARY_PATH="${CRAY_MPICH_DIR}/lib:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${CRAY_MPICH_DIR}/lib:${LD_LIBRARY_PATH}"

export CXX=hipcc
export CXXFLAGS+=" -ferror-limit=1"

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=7,1

export ESPRESO_RANK_TO_GPU_MAP="0"
# export ESPRESO_RANK_TO_GPU_MAP="4,5,2,3,6,7,0,1" # dont use

export ESPRESO_ROCM_ARCH="gfx90a:sramecc+:xnack-"

export ESPRESO_USE_WRAPPER_DNBLAS=mkl
export ESPRESO_USE_WRAPPER_DNSOLVER=mkl
export ESPRESO_USE_WRAPPER_LAPACK=mkl
export ESPRESO_USE_WRAPPER_SPBLAS=mkl
export ESPRESO_USE_WRAPPER_SPSOLVER=mkl
export ESPRESO_USE_WRAPPER_SCSOLVER=mkl
export ESPRESO_USE_WRAPPER_GPU=rocm



# always use slurm options --gpus-per-task=1 --cpus-per-task=7
# salloc --account=project_465000572 --partition=standard-g --ntasks=8 --gpus-per-task=1 --cpus-per-task=7  --time=8:00:00
# srun -n 8 ./build/espreso -c path/to/espreso.ecf

#!/usr/bin/env bash



# build-institution-machine-compiler-mpi-blas-specialname
buildname="build-anl-aurora-intel-mpich-mkl-oneapi"
export WAFLOCK=".lock-waf_linux_${buildname}"
curr_buildname="$(readlink build)"
if [ "${curr_buildname}" != "${buildname}" ]
then
    rm -rf build
    ln -s "${buildname}" build
fi



module purge

module load oneapi/release/2025.0.5
module load mpich/opt/develop-git.6037a7a
module load cray-pals/1.4.0

if [ ! -d "dependencies/oneAPI-2025.2.1" ]
then
    echo "intel toolkit not installed in dependencies, please install it manually. DO NOT install mpi"
    return 3
fi

setvars_config=$(mktemp)
echo "mpi=exclude" > "${setvars_config}"
source dependencies/oneAPI-2025.2.1/setvars.sh --config="${setvars_config}"
rm -f "${setvars_config}"

export CPATH="${PWD}/dependencies/oneAPI-2025.2.1/2025.2/include:${CPATH}"

export MPICH_CC=icx
export MPICH_CXX=icpx
export MPICH_FC=ifx



. env/dependencies/install.cmake.sh intel_mpich_mkl_oneapi x86_64
. env/dependencies/install.elfutils.sh intel_mpich_mkl_oneapi icx
. env/dependencies/install.gklib.sh intel_mpich_mkl_oneapi icx
. env/dependencies/install.metis32.sh intel_mpich_mkl_oneapi icx
. env/dependencies/install.parmetis32.sh intel_mpich_mkl_oneapi mpicc
. env/dependencies/install.suitesparse.sh intel_mpich_mkl_oneapi icx ifx
. env/dependencies/install.pastix.sh intel_mpich_mkl_oneapi icpx icx ifx Intel10_64lp_seq



export CXX=mpicxx
export ES_INT_WIDTH=32
export CXXFLAGS+=" -fmax-errors=1"
export CXXFLAGS+=" -DESPRESO_STACKTIMER_ENABLE"

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=8,1
export CPU_BIND_SCHEME_12ppn="--cpu-bind=list:1-8:9-16:17-24:25-32:33-40:41-48:53-60:61-68:69-76:77-84:85-92:93-100"

# export ESPRESO_MAP_LOCALRANK_TO_GPU="$(seq -s' ' 0 11)"

export ZE_FLAT_DEVICE_HIERARCHY="FLAT" # two stacks on a GPU card appear each as a single GPU device
export ZES_ENABLE_SYSMAN="1" # to make free_memory device info available, enabled by default since 2025.something

export ESPRESO_SYCL_TARGETS="spir64_gen"
export ESPRESO_SYCL_BACKEND_OPTIONS="-device pvc"

export ESPRESO_USE_WRAPPER_DNBLAS=mkl
export ESPRESO_USE_WRAPPER_DNSOLVER=mkl
export ESPRESO_USE_WRAPPER_LAPACK=mkl
export ESPRESO_USE_WRAPPER_SPBLAS=mkl
export ESPRESO_USE_WRAPPER_SPSOLVER=mkl
export ESPRESO_USE_WRAPPER_SCSOLVER=mkl
export ESPRESO_USE_WRAPPER_GPU=oneapi



# mpirun ${CPU_BIND_SCHEME_12ppn} -ppn 12 -n XX ./app

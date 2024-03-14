#!/bin/bash

DEPENDENCIES_DIR="dependencies"
mkdir -p dependencies

module use /pfs/lustrep2/projappl/project_462000125/samantao-public/mymodules
module load LUMI/23.09
module load partition/G
module load rocm/5.4.3



AOCL_ROOT="${DEPENDENCIES_DIR}/aocl-linux-aocc-4.1.0"
if [ ! -d "${AOCL_ROOT}" ]
then
    (
        echo "AOCL not installed, installing locally..."
        cd "${DEPENDENCIES_DIR}"
        curl -s -O https://download.amd.com/developer/eula/aocl/aocl-4-1/aocl-linux-aocc-4.1.0.tar.gz
        tar -xf aocl-linux-aocc-4.1.0.tar.gz
        rm aocl-linux-aocc-4.1.0.tar.gz
        cd aocl-linux-aocc-4.1.0
        ./install.sh -t "${PWD}" -i lp64
    )
fi
source "${AOCL_ROOT}/4.1.0/aocc/amd-libs.cfg"



METIS_ROOT="${DEPENDENCIES_DIR}/metis-5.1.0"
if [ ! -d "${METIS_ROOT}" ]
then
    (
        echo "METIS not installed, installing locally..."
        cd "${DEPENDENCIES_DIR}"
        curl -s -O http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
        tar -xf metis-5.1.0.tar.gz
        rm metis-5.1.0.tar.gz
        cd metis-5.1.0
        sed -i 's/#define REALTYPEWIDTH 32/#define REALTYPEWIDTH 64/' include/metis.h
        make config shared=1 cc=clang gklib_path="${PWD}/GKlib" prefix="${PWD}"
        make -j$(nproc)
        make install
    )
fi
export CPATH="${PWD}/${METIS_ROOT}/include:${CPATH}"
export LIBRARY_PATH="${PWD}/${METIS_ROOT}/lib:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH="${PWD}/${METIS_ROOT}/lib:${LD_LIBRARY_PATH}"



VERSION_SUITESPARSE=v7.6.0
export SUITESPARSE_DIR="SuiteSparse_${VERSION_SUITESPARSE}_clang_aocl"
export SUITESPARSE_ROOT="${DEPENDENCIES_DIR}/${SUITESPARSE_DIR}"
if [ ! -d "${SUITESPARSE_ROOT}" ]
then
    (
        echo "SuiteSparse ${VERSION_SUITESPARSE} not installed, installing locally..."
        cd "${DEPENDENCIES_DIR}"
        git clone -b "${VERSION_SUITESPARSE}" https://github.com/DrTimothyAldenDavis/SuiteSparse.git "${SUITESPARSE_DIR}"
        cd "${SUITESPARSE_DIR}"

        module load buildtools/23.09
        mkdir -p build
        cd build
        cmake -DCMAKE_C_COMPILER=clang -DENABLE_CUDA=false -DSUITESPARSE_ENABLE_PROJECTS="cholmod" -DLAPACK_LIBRARIES="${DEPENDENCIES_DIR}/aocl-linux-aocc-4.1.0/4.1.0/aocc/lib_LP64/libflame.so.4.1.0" -DCMAKE_INSTALL_PREFIX=.. ..
        cmake --build . -j $(nproc)
        cmake --install .
    )
fi
export CPATH="${PWD}/${SUITESPARSE_ROOT}/include:${CPATH}"
export LIBRARY_PATH="${PWD}/${SUITESPARSE_ROOT}/lib64:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH="${PWD}/${SUITESPARSE_ROOT}/lib64:${LD_LIBRARY_PATH}"

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



# always use slurm options --gpus-per-task=1 --cpus-per-task=7
# salloc --account=project_465000572 --partition=standard-g --ntasks=8 --gpus-per-task=1 --cpus-per-task=7  --time=8:00:00
# srun -n 8 ./build/espreso -c path/to/espreso.ecf

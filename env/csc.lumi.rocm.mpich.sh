#!/bin/bash



module load LUMI/24.03
module load partition/G
module load wget/1.21.4-cpeGNU-24.03
module load rocm/6.0.3
module load aocc/4.1.0



. dependencies/install.cmake.sh
. dependencies/install.aocl.sh
. dependencies/install.suitesparse.sh hipcc hipfc "-DLAPACK_LIBRARIES=${PWD}/dependencies/aocl-linux-aocc-4.1.0/install/4.1.0/aocc/lib/libflame.so.4.1.0"
. dependencies/install.gklib.sh hipcc
. dependencies/install.metis32.sh hipcc
. dependencies/install.parmetis32.sh mpicc



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



# always use slurm options --gpus-per-task=1 --cpus-per-task=7
# salloc --account=project_465000572 --partition=standard-g --ntasks=8 --gpus-per-task=1 --cpus-per-task=7  --time=8:00:00
# srun -n 8 ./build/espreso -c path/to/espreso.ecf

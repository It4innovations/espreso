
ml Python/3.10.4-GCCcore-11.3.0
ml OpenMPI/4.1.4-GCC-11.3.0
ml OpenBLAS/0.3.20-GCC-11.3.0
ml CMake/3.23.1-GCCcore-11.3.0

. dependencies/install.suitesparse.sh gcc
. dependencies/install.metis32.sh gcc
. dependencies/install.parmetis32.sh mpicc

export CXX=mpic++
export BLAS_LIBRARIES=openblas
export LAPACK_LIBRARIES=openblas
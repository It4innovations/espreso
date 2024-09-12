
ml OpenMPI/4.1.6-GCC-13.2.0
ml OpenBLAS/0.3.24-GCC-13.2.0
ml imkl/2024.2.0
ml CMake/3.27.6-GCCcore-13.2.0
ml Boost/1.83.0-GCC-13.2.0
ml Eigen/3.4.0-GCCcore-13.2.0
ml Python/3.11.5-GCCcore-13.2.0
ml Boost.Python-NumPy

. dependencies/install.suitesparse.sh gcc
. dependencies/install.gklib.sh gcc
. dependencies/install.metis32.sh gcc
. dependencies/install.parmetis32.sh mpicc
. dependencies/install.precice.sh mpic++

export CXX=mpic++
export BLAS_LIBRARIES=openblas
export LAPACK_LIBRARIES=openblas

export OMPI_MCA_btl=^openib
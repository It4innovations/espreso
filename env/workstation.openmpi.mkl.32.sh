
buildname="build-workstation-openmpi-mkl-32"
export WAFLOCK=".lock-waf_linux_${buildname}"
rm -rf build
ln -s "${buildname}" build



. env/dependencies/install.suitesparse.sh mkl gcc gfortran
. env/dependencies/install.gklib.sh mkl gcc
. env/dependencies/install.metis32.sh mkl gcc
. env/dependencies/install.parmetis32.sh mkl mpicc
. env/dependencies/install.precice.sh mkl mpic++
. env/dependencies/install.mumps.sh mkl mpicc mpifort "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core" "-lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core" "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core"
. env/dependencies/install.strumpack.sh mkl g++ gcc gfortran "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core" "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core"
. env/dependencies/install.pastix.sh mkl g++ gcc
. env/dependencies/install.superlu_dist.sh mkl g++ gcc "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core"

export CXX=mpic++
export ES_INT_WIDTH=32
export CPATH=/usr/include/mkl:$CPATH

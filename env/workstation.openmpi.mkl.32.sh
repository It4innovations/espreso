
. env/dependencies/install.suitesparse.sh mkl gcc gfortran
. env/dependencies/install.gklib.sh mkl gcc
. env/dependencies/install.metis32.sh mkl gcc
. env/dependencies/install.parmetis32.sh mkl mpicc
. env/dependencies/install.precice.sh mkl mpic++

export CXX=mpic++
export ES_INT_WIDTH=32
export CPATH=/usr/include/mkl:$CPATH

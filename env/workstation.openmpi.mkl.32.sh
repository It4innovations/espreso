
. env/dependencies/install.suitesparse.sh gcc gfortran
. env/dependencies/install.gklib.sh gcc
. env/dependencies/install.metis32.sh gcc
. env/dependencies/install.parmetis32.sh mpicc
. env/dependencies/install.precice.sh mpic++

export CXX=mpic++
export ES_INT_WIDTH=32
export CPATH=/usr/include/mkl:$CPATH

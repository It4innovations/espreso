
. dependencies/install.suitesparse.sh gcc
. dependencies/install.gklib.sh gcc
. dependencies/install.metis32.sh gcc
. dependencies/install.parmetis32.sh mpicc
. dependencies/install.precice.sh mpic++

export CXX=mpic++
export ES_INT_WIDTH=32
export CPATH=/usr/include/mkl:$CPATH

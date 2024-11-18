
. dependencies/install.suitesparse.sh gcc
. dependencies/install.gklib.sh gcc
. dependencies/install.metis64.sh gcc
. dependencies/install.parmetis64.sh mpicc
. dependencies/install.precice.sh mpic++

export CXX=mpic++
export CPATH=/usr/include/mkl:$CPATH

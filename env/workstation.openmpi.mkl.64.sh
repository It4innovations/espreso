
. dependencies/install.metis64.sh gcc
. dependencies/install.parmetis64.sh mpicc

export CXX=mpic++
export CPATH=/usr/include/mkl:$CPATH
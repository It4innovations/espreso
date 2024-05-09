
. dependencies/install.metis32.sh gcc
. dependencies/install.parmetis32.sh mpicc

export CXX=mpic++
export CPATH=/usr/include/mkl:$CPATH

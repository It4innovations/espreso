
. dependencies/install.gklib.sh gcc
. dependencies/install.metis32.sh gcc
. dependencies/install.parmetis32.sh mpicc
. dependencies/install.precice.sh mpic++

export CXX=mpic++
export CPATH=/usr/include/mkl:$CPATH

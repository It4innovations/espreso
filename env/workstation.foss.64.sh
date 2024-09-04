
. dependencies/install.suitesparse.sh gcc gfortran
. dependencies/install.gklib.sh gcc
. dependencies/install.metis64.sh gcc
. dependencies/install.parmetis64.sh mpicc
. dependencies/install.precice.sh mpic++

export CXX=mpic++
export ES_INT_WIDTH=64
